! ------------------------------------------------------------------------------
! ATS Driver module
!
! Effectively a singleton module for ATS points of entry.
!
! Author: Ethan Coon (ecoon _at_ lanl.gov)
! ------------------------------------------------------------------------------

module ats_clm
  use clm_precision
  use clm_type_module
  use clm_host
  use clm_host_transfer
  use clm1d_type_module
  use io_type_module
  implicit none

  private

  ! namespace-local data to hold the singleton clm and host instances  
  type(clm_type) :: clm
  type(host_type) :: host
  real(r8) :: year0


  public :: ats_to_clm_ground_properties, &
       ats_to_clm_dz, &
       ats_to_clm_met_data, &
       ats_to_clm_forced_vegetation, &
       ats_to_clm_pressure, &
       ats_to_clm_wc, &
       ats_to_clm_tksat_from_porosity, &
       ats_to_clm_irrigation, &
       ats_to_clm_et_controls, &
       clm_to_ats_ground_energy_fluxes, &
       clm_to_ats_total_energy_fluxes, &
       clm_to_ats_mass_fluxes, &
       clm_to_ats_total_mass_fluxes, &
       clm_to_ats_total_mass_fluxes_combined, &
       clm_to_ats_diagnostics, &
       ats_clm_init, &
       ats_clm_setup_begin, &
       ats_clm_setup_end, &
       ats_clm_advance_time

contains


  !
  ! Begin initialization, allocating space for driver, grid.
  ! ------------------------------------------------------------------
  ! Input:
  !   ncells     | number of grid cells in the subsurface
  !   ncolumns   | number of columns.  Note it is currently assumed that
  !              |  ncolumns divides ncells evenly, with a fixed number of
  !              |  cells per column.  (FIXME)
  !   startcode  | 1 if new run, 2 if restart (not currently supported?)
  !   rank       | rank (for logging purposes only)
  !   verbosity  | 0 (None) - 1 (Low) - 2 (High) - 3 (Extreme)
  !   
  subroutine ats_clm_init(ncells, ncolumns, col_inds, startcode, rank, verbosity) bind(C)
    use clm_io_config,only : MAX_FILENAME_LENGTH
    implicit none
    integer(i4),intent(in) :: ncells, ncolumns, startcode, rank, verbosity
    integer(i4),intent(in) :: col_inds(ncolumns, 2)

    ! locals
    character(len=MAX_FILENAME_LENGTH) :: output_dir = "clm/"

    !--- initialize the clm master object
    call clm_init(clm, rank, ncolumns, 1, 18, verbosity)

    !--- open log files, set up io options
    print*, "Opening logfile: ", output_dir
    call io_open(clm%io, output_dir, rank, 0)
    clm%io%restart_last = 1
    clm%io%restart_daily = 0
    clm%io%dump_interval = -1
    clm%io%dump_current = 0
    clm%io%output_1d = 0
    clm%io%write_binaries = 0

    !--- initialize the host object, pushing grid info into it
    call host_init(host, ncells, ncells, ncolumns, ncolumns, col_inds)
    if (io_ok(clm%io, VERBOSITY_LOW)) call host_write_to_log(host, clm%io%log)

    !--- read clm input files, data
    ! FIXME -- can we remove this and use a setter? --etc
    ! FIXME -- fix io to not pass clm_write_logs, instead iounit --etc
    clm%drv%maxt = 1  ! hard-coded 1 PFT per cell (FIXME)
    clm%drv%vclass = 2 ! likely can be removed completely (FIXME)
    clm%drv%startcode = startcode
    clm%drv%clm_ic = startcode

    if (io_ok(clm%io, VERBOSITY_LOW)) then
       write(clm%io%log,*) "  CLM startcode for date (1=restart, 2=defined):", clm%drv%startcode
       write(clm%io%log,*) "  CLM IC (1=restart, 2=defined):", clm%drv%clm_ic
       if (clm%drv%startcode == 0) then
          write(clm%io%log,*) "ERROR: startcode = 0"
          stop
       end if
       if (clm%drv%clm_ic == 0) then
          write(clm%io%log,*) "ERROR: clm_ic = 0"
          stop
       end if
    end if

    clm%clm => clm1d_create_n(clm%ntiles, clm%drv%surfind, clm%drv%soilind, clm%drv%snowind)
  end subroutine ats_clm_init


  !
  ! Set the reference, time 0, in years.
  ! ------------------------------------------------------------------
  ! Input:
  !   zero_time  | The time from which all time is measured. [year]
  !              |  i.e. 2010.0 for midnight Jan 1, 2010.  Used in
  !              |  determining sun/shade factors for canopy.
  !   
  subroutine ats_clm_zero_time(zero_year) bind(C)
    ! zero time, in years
    implicit none
    real(r8),intent(in) :: zero_year
    year0 = zero_year
    call drv_time2date(year0, clm%drv%doy, clm%drv%day, clm%drv%gmt,       &
         clm%drv%yr, clm%drv%mo, clm%drv%da, clm%drv%hr,        &
         clm%drv%mn, clm%drv%ss)
    call drv_time2date(year0, clm%drv%sdoy, clm%drv%day, clm%drv%sgmt,       &
         clm%drv%syr, clm%drv%smo, clm%drv%sda, clm%drv%shr,        &
         clm%drv%smn, clm%drv%sss)
  end subroutine ats_clm_zero_time


  !
  ! Sets the initial, assumed uniform, state.
  ! ------------------------------------------------------------------
  ! Input:
  !   temperature| Soil, snow, and water temperature [K]
  !   snow depth | Initial snow depth [m]
  !   
  subroutine ats_clm_initial_state(temperature, snow_depth) bind(C)
    ! K, [mm]
    implicit none
    real(r8),intent(in) :: temperature  ! uniform initial temp [K]
    real(r8),intent(in) :: snow_depth   ! uniform initial snow depth [m]
    clm%drv%t_ini = temperature
    clm%drv%h2osno_ini = snow_depth * 0.1 * 1000. ! CLM expects SWE [mm]
  end subroutine ats_clm_initial_state


  !
  ! Set soil properties
  ! ------------------------------------------------------------------
  ! Input:
  !   latlon     | lat/lon.  +Northern and Eastern hemispheres. [degrees]
  !              |  Array of shape [ncolumns, 2]
  !   sand,clay  | soil texture, fractions (must range from 0-1) [-] Size ncells.
  !   color_index| index into soil color models?  See clm1d_varcon albsat and albdry
  !              |  Size ncolumns.
  !   f_ground   | Fraction of each land type.  Land types are set in drv_vegp.dat
  !              |  Array of shape [ncolumns, NUM_LC_CLASSES]
  !              |  Note only the dominant land type is used, as maxt is
  !              |  hard-coded 1 and this assumption is used in a few places in drv.
  !   
  subroutine ats_to_clm_ground_properties(latlon, sand, clay, &
       color_index, fractional_ground) bind(C)
    implicit none
    real(r8),intent(in) :: latlon(2,host%ncolumns_g)    ! latitude,longitude [decimal degrees]
    real(r8),intent(in) :: sand(host%ncells_g)          ! percent sand FIXME: 0-1 or 0-100? --etc
    real(r8),intent(in) :: clay(host%ncells_g)          ! percent clay FIXME: 0-1 or 0-100? --etc
    integer(i4),intent(in) :: color_index(host%ncolumns_g)  ! color index FIXME: document! --etc
    real(r8),intent(in) :: fractional_ground(clm%drv%nt,host%ncolumns_g) ! fraction of land surface of type t
    call host_to_clm_ground_properties(host, latlon, sand, clay, color_index, fractional_ground, clm)
  end subroutine ats_to_clm_ground_properties


  !
  ! Begin setup. Sets the tiles and pushes info from driver into grid/tile
  ! ------------------------------------------------------------------
  !   
  subroutine ats_clm_setup_begin() bind(C)
    implicit none
    call clm_setup_begin(clm)
  end subroutine ats_clm_setup_begin


  !
  ! sets the cell thicknesses in the vertical
  ! ------------------------------------------------------------------
  ! Input:
  !   dz         | vector of cell dz [m]
  !
  subroutine ats_to_clm_dz(dz) bind(C)
    implicit none
    real(r8),intent(in) :: dz(host%ncells_g)
    call host_to_clm_dz(host, 1.d0, dz, clm)
    clm%clm(:)%soi_z = 5 ! FIXME made up, refernce temperature cell.  should be function of dz
  end subroutine ats_to_clm_dz


  !
  ! Sets ET controls
  ! ------------------------------------------------------------------
  ! Input:
  !   beta_type          | ??
  !   veg_water_stress_type
  !                      | ??
  !   wilting_point      | ??
  !   field_capacity     | ??
  !   res_sat            | residual saturation [-]
  !                      |  make this variable across columns?  cells? FIXME
  !              
  subroutine ats_to_clm_et_controls(beta_type, veg_water_stress_type, wilting_point, &
       field_capacity, res_sat) bind(C)
    implicit none

    ! ET controls
    integer(i4), intent(in) :: beta_type              ! beta formulation for bare soil Evap 0=none, 1=linear, 2=cos
    integer(i4), intent(in) :: veg_water_stress_type  ! veg transpiration water stress formulation
    ! 0=none, 1=press, 2=sm
    real(r8), intent(in) :: wilting_point         ! wilting point in m if press-type, in saturation
    ! if soil moisture type
    real(r8), intent(in) :: field_capacity        ! field capacity for water stress same as units above
    real(r8), intent(in) :: res_sat               ! residual saturation
    call host_to_clm_et_controls(host, beta_type, veg_water_stress_type, wilting_point, &
         field_capacity, res_sat, clm)
  end subroutine ats_to_clm_et_controls


  !
  ! End setup. Pushes grid, tile, drv info into clm1d column instances.
  ! ------------------------------------------------------------------
  !   
  subroutine ats_clm_setup_end() bind(C)
    implicit none
    call clm_setup_end(clm)
  end subroutine ats_clm_setup_end


  !
  ! Move water content data into clm1d columns
  ! ------------------------------------------------------------------
  ! Input:
  !   porosity   | [-]  Size ncells.
  !   saturation | [-]  Size ncells.
  !
  subroutine ats_to_clm_wc(porosity, saturation) bind(C)
    implicit none
    real(r8),intent(in) :: porosity(host%ncells_g)      ! [-]
    real(r8),intent(in) :: saturation(host%ncells_g)    ! [-]
    call host_to_clm_wc(host, porosity, saturation, clm)
  end subroutine ats_to_clm_wc


  !
  ! Set saturated thermal conductivity from porosity.
  ! -----------------------------------------------------------------------
  !  NOTE: This is seperate from set_wc() to allow it to be called once and
  !   assume to be fixed throughout the simulation (small porosity
  !   variation due to compressibility only.
  !
  ! Input:
  !   porosity   | [-]  Size ncells.  
  !
  subroutine ats_to_clm_tksat_from_porosity(poro) bind(C)
    implicit none
    real(r8),intent(in) :: poro(host%ncells_g)  ! [-]
    call host_to_clm_tksat_from_porosity(host, poro, clm)
  end subroutine ats_to_clm_tksat_from_porosity


  !
  ! Move pressure data into 1d columns.
  ! ------------------------------------------------------------------
  ! Input:
  !   pressure   | [Pa] Size ncells.
  !   p_atm      | Atmospheric pressure.  [Pa] Note CLM works in head
  !              |  units, so p_atm is used to convert Pa to mm.  It
  !              |  isn't obvious that this is a good idea in cases of
  !              |  frozen water?  FIXME
  !
  subroutine ats_to_clm_pressure(pressure, p_atm) bind(C)
    use clm1d_varcon, only : denh2o, grav
    implicit none
    real(r8),intent(in) :: pressure(host%ncells_g)  ! [Pa]
    real(r8),intent(in) :: p_atm

    ! local
    real(r8) :: pressure_adj(host%ncells_g)  ! [mm] 
    pressure_adj(:) = (pressure(:) - p_atm) / (denh2o*grav) * 1e3
    call host_to_clm_pressure(host, pressure_adj, 1.d0, clm)
  end subroutine ats_to_clm_pressure


  !
  ! Set the meteorological data.
  ! ------------------------------------------------------------------
  ! NOTE: all are of size ncolumns.
  !
  ! NOTE: this is a little wonky, as CLM1D actually wants split precip,
  !  and ATS uses split precip, but the CLM driver expects summed
  !  precip that it splits by air temperature.  FIXME
  !
  ! Input:
  !   qSW        | Shortwave incoming radiation [W/m^2] Size ncolumns.
  !   qLW        | Longwave incoming radiation [W/m^2] Size ncolumns.
  !   pRain      | Rainfall precipitation rate [m/s]
  !   pSnow      | Snowfall precipitation rate [m/s]
  !   air_temp   | Air temperature [K]
  !   vp_air     | Vapor presure air [Pa]
  !   wind_u     | Windspeed velocity [m/s]
  !   p_atm      | Atmospheric pressure [Pa]
  !
  subroutine ats_to_clm_met_data(eflx_swin, eflx_lwin, precip, &
       air_temp, vp_air, wind_x, wind_y, patm) bind(C)
    implicit none
    real(r8),intent(in) :: eflx_swin(host%ncolumns_g) ! shortwave incoming radiation [W/m^2]
    real(r8),intent(in) :: eflx_lwin(host%ncolumns_g) ! longwave incoming radiation [W/m^2]
    real(r8),intent(in) :: precip(host%ncolumns_g) ! precipitation rate [mm/s]
    real(r8),intent(in) :: air_temp(host%ncolumns_g) ! air temperature [K]
    real(r8),intent(in) :: vp_air(host%ncolumns_g) ! vapor pressure air [Pa]
    real(r8),intent(in) :: wind_x(host%ncolumns_g) ! wind speed, eastward direction [m/s]
    real(r8),intent(in) :: wind_y(host%ncolumns_g) ! wind speed, northward direction [m/s]
    real(r8),intent(in) :: patm(host%ncolumns_g) ! atmospheric pressure [Pa]
    
    ! local conversion to specific humidity [kg/kg]
    real(r8) :: spec_hum(host%ncolumns_g)
    integer i
    
    do i=1,host%ncolumns_g
       ! convert to specific humidity
       spec_hum(i) = 0.622 * vp_air(i) / (patm(i) - 0.378 * vp_air(i))
    end do

    call host_to_clm_met_data(host, eflx_swin, eflx_lwin, precip, &
         air_temp, spec_hum, wind_x, wind_y, patm, clm)
  end subroutine ats_to_clm_met_data



  !
  ! Sets LAI, SAI for Sattelite/fixed phenology mode.
  ! ------------------------------------------------------------------
  subroutine ats_to_clm_forced_vegetation(lai, sai, z0m, displacement_ht)
    implicit none
    real(r8),intent(in) :: lai(host%ncolumns_g) ! exposed leaf area index [-]
    real(r8),intent(in) :: sai(host%ncolumns_g) ! exposed stem area index [-]
    real(r8),intent(in) :: z0m(host%ncolumns_g) ! aerodynamic roughness length [m]
    real(r8),intent(in) :: displacement_ht(host%ncolumns_g) ! displacement height [m]
    call host_to_clm_forced_vegetation(host, lai, sai, z0m, displacement_ht, clm)
  end subroutine ats_to_clm_forced_vegetation


  !
  ! Sets irrigiation data
  !
  ! Expected units: 
  ! ------------------------------------------------------------------
  subroutine ats_to_clm_irrigation(irr_type, irr_cycle, &
       irr_rate, irr_start, irr_stop, irr_threshold, irr_thresholdtype) bind(C)
    implicit none
    ! irrigation keys
    integer(i4), intent(in) :: irr_type            ! irrigation type flag (0=none,1=spray,2=drip,3=instant)
    integer(i4), intent(in) :: irr_cycle           ! irrigation cycle flag (0=constant,1=deficit)
    real(r8), intent(in) :: irr_rate           ! irrigation application rate for spray and drip [mm/s]
    real(r8), intent(in) :: irr_start          ! irrigation daily start time for constant cycle
    real(r8), intent(in) :: irr_stop           ! irrigation daily stop tie for constant cycle
    real(r8), intent(in) :: irr_threshold      ! irrigation threshold criteria for deficit cycle
    ! (units of soil moisture content)
    integer(i4), intent(in)  :: irr_thresholdtype  ! irrigation threshold criteria type -- top layer,
    ! bottom layer, column avg
    call host_to_clm_irrigation(host, irr_type, irr_cycle, &
         irr_rate, irr_start, irr_stop, irr_threshold, irr_thresholdtype, clm)
  end subroutine ats_to_clm_irrigation


  !
  ! Advance the timestep
  ! ------------------------------------------------------------------
  ! Input:
  !   step       | integer cycle number (logging only?)
  !   time       | time at start of step (relative to zero time) [s]
  !   dt         | step size [s]
  subroutine ats_clm_advance_time(istep, time, dtime) bind(C)
    implicit none
    integer(i4),intent(in) :: istep
    real(r8),intent(in) :: time, dtime
    if ((istep.eq.0).or.(clm%drv%time < 0)) then
       ! tick to set initial time
       clm%drv%ts = nint(time)
       call drv_tick(clm%drv)
    end if
    if (io_ok(clm%io, VERBOSITY_LOW)) then
       write(clm%io%log,*) "Advancing: step ", istep, " from time ", time, " with size ", dtime
    end if
    call clm_advance_time(clm, host, istep, time, dtime)

    print*, "DynVegPar: elai = ", clm%clm%tlai, clm%clm%elai
  end subroutine ats_clm_advance_time


  !
  ! Source/sink terms for integrated hydrology code.
  ! ------------------------------------------------------------------
  ! Output:
  !   qW_surf    | water source/sink surface (sign?) [mm/s] Size ncolumns.
  !   qW_subsurf | water source/sink surface (sign?) [mm/m/s] Size ncells.
  !
  subroutine clm_to_ats_total_mass_fluxes(qflx_surface, qflx_subsurface) bind(C)
    implicit none
    real(r8),intent(out) :: qflx_surface(host%ncolumns_g)  ! total mass flux to/from the surface [mm/s]
    real(r8),intent(out) :: qflx_subsurface(host%ncells_g) ! total mass flux to/from the subsurface [1/s]
    call clm_to_host_total_mass_fluxes(host, clm, qflx_surface, qflx_subsurface)
  end subroutine clm_to_ats_total_mass_fluxes



  !
  ! Energy fluxes for diagnostics/visualization
  ! ------------------------------------------------------------------
  ! Note: all sizes are ncolumns
  !
  ! Output:
  !   latent_heat        | Latent heat flux [W/m^2] (sign?)
  !   sensible_heat      | Sensible heat flux [W/m^2] (sign?)
  !   longwave_out       | Outward longwave radiation from surface to atmosphere [W/m^2]
  !   conducted_e        | Energy conducted to subsurface (sign?) [W/m^2]
  !
  subroutine clm_to_ats_total_energy_fluxes(eflx_lh, eflx_sh, eflx_lwrad_out, eflx_soil) bind(C)
    implicit none
    real(r8),intent(out) :: eflx_lh(host%ncolumns_g)          ! latent heat flux [W/m^2]
    real(r8),intent(out) :: eflx_sh(host%ncolumns_g)          ! sensible heat from ground [W/m^2]
    real(r8),intent(out) :: eflx_lwrad_out(host%ncolumns_g)   ! outgoing long-wave radiation from ground [W/m^2]
    real(r8),intent(out) :: eflx_soil(host%ncolumns_g)        ! flux conducted to ground [W/m^2]
    call clm_to_host_total_energy_fluxes(host, clm, eflx_lh, eflx_sh, eflx_lwrad_out, eflx_soil)
    eflx_sh(:) = -eflx_sh(:) ! sensible + to ground + leaf
    eflx_lh(:) = -eflx_lh(:) ! latent + to ground + leaf
  end subroutine clm_to_ats_total_energy_fluxes


  !
  ! Ground components (not canopy) of energy fluxes
  ! ------------------------------------------------------------------
  ! Note: all sizes are ncolumns
  !
  ! Output:
  !   latent_heat        | Latent heat flux [W/m^2] (+ to leaf/ground)
  !   sensible_heat      | Sensible heat flux [W/m^2] (+ to leaf/ground)
  !   longwave_out       | longwave radiation [W/m^2] (+ to atmosphere)
  !   conducted_e        | Energy conducted to subsurface [W/m^2] (+ to soil)
  !
  subroutine clm_to_ats_ground_energy_fluxes(eflx_lh, eflx_sh, eflx_lwrad_out, eflx_soil) bind(C)
    implicit none
    real(r8),intent(out) :: eflx_lh(host%ncolumns_g)          ! latent heat flux [W/m^2]
    real(r8),intent(out) :: eflx_sh(host%ncolumns_g)          ! sensible heat from ground [W/m^2]
    real(r8),intent(out) :: eflx_lwrad_out(host%ncolumns_g)   ! outgoing long-wave radiation from ground [W/m^2]
    real(r8),intent(out) :: eflx_soil(host%ncolumns_g)        ! flux conducted to ground [W/m^2]
    call clm_to_host_ground_energy_fluxes(host, clm, eflx_lh, eflx_sh, eflx_lwrad_out, eflx_soil)
    eflx_sh(:) = -eflx_sh(:) ! sensible + to ground
    eflx_lh(:) = -eflx_lh(:) ! latent + to ground
  end subroutine clm_to_ats_ground_energy_fluxes



  !
  ! Mass fluxes for mass balance and diagnostics/visualization
  ! ------------------------------------------------------------------
  ! Note: size ncolumns except where otherwise noted.
  !
  ! Output:
  !   evap_total         | Total evaporation to atmosphere [mm/s]
  !   evap_ground        | Ground component of evaporation (no snow sublimation) [mm/s]
  !   evap_soil          | Soil evaporation [mm/s]
  !   evap_canopy        | Canopy component of evaporation [mm/s]
  !   tran_veg           | Transpiration from vegetation over the column [mm/s]
  !   influx             | Precip? ???? [mm/s]
  !   irrigation         | Surface irrigation (can be intercepted) [mm/s]
  !   irrigation_inst    | Irrigation delivered to subsurface directly [mm/m/s???] (Size ncells)
  !   irrigation_flag    | ???
  !   tran_soil          | Transpiration distributed to soil via rooting curve [mm/m/s] (Size ncells)
  !
  subroutine clm_to_ats_mass_fluxes(qflx_evap_tot, qflx_evap_ground, qflx_evap_soil, &
       qflx_evap_veg, qflx_tran_veg, qflx_infl, qflx_irr, qflx_irr_inst, irr_flag, &
       qflx_tran_soil) bind(C)
    implicit none
    real(r8),intent(out) :: qflx_evap_tot(host%ncolumns_g)    ! total evaporation [mm/s]
    real(r8),intent(out) :: qflx_evap_ground(host%ncolumns_g) ! ground evaporation (does not
    !  include snow sublimation) [mm/s]
    real(r8),intent(out) :: qflx_evap_soil(host%ncolumns_g)   ! soil evaporation [mm/s]
    real(r8),intent(out) :: qflx_evap_veg(host%ncolumns_g)    ! canopy evaporation [mm/s]
    real(r8),intent(out) :: qflx_tran_veg(host%ncolumns_g)    ! canopy transpiration [mm/s]
    real(r8),intent(out) :: qflx_infl(host%ncolumns_g)        ! net infiltration [mm/s]
    real(r8),intent(out) :: qflx_irr(host%ncolumns_g)         ! irrigation sources [mm/s]
    real(r8),intent(out) :: qflx_irr_inst(host%ncells_g)      ! instantaneous irrigation sources [mm/s]
    real(r8),intent(out) :: irr_flag(host%ncolumns_g)         ! flag for irrigation type
    real(r8),intent(out) :: qflx_tran_soil(host%ncells_g)     ! transpiration, distributed via roots [mm/s]
    call clm_to_host_mass_fluxes(host, clm, qflx_evap_tot, qflx_evap_ground, qflx_evap_soil, &
         qflx_evap_veg, qflx_tran_veg, qflx_infl, qflx_irr, qflx_irr_inst, irr_flag, qflx_tran_soil)
  end subroutine clm_to_ats_mass_fluxes



  !
  ! Source/sink terms for integrated hydrology code, with surface fluxes
  ! added into the top cell of the subsurface.
  ! ------------------------------------------------------------------
  ! Output:
  !   qW_subsurf | water source/sink surface (sign?) [mm/m/s] Size ncells.
  !
  subroutine clm_to_ats_total_mass_fluxes_combined(qflx_subsurface) bind(C)
    implicit none
    real(r8),intent(out) :: qflx_subsurface(host%ncells_g) ! total mass flux to/from the subsurface [1/s]
    call clm_to_host_total_mass_fluxes_combined(host, clm, qflx_subsurface)
  end subroutine clm_to_ats_total_mass_fluxes_combined


  !
  ! Diagnostic variables
  ! ------------------------------------------------------------------
  ! Note: size ncolumns except where otherwise noted.
  !
  ! Output:
  !   swe                | snow water equivalent [m]
  !   snow_depth         | [m] CHECK THIS!  CLM thinks this is [m] but everything
  !                      |    else is [mm]... FIXME
  !   canopy_storage     | water stored in the canopy [m]
  !   T_skin             | Surface skin temperature [K]
  !   T_veg              | Leaf temperature [K]
  !   T_soil             | Soil temperature [K]  (size ncells)
  !
  subroutine clm_to_ats_diagnostics(swe, snow_depth, canopy_storage, T_skin, T_veg, T_soil) bind(C)
    implicit none
    real(r8),intent(out) :: swe(host%ncolumns_g) ! snow-water equivalent (mass) [kg/m^2]
    real(r8),intent(out) :: canopy_storage(host%ncolumns_g) ! canopy storage (mass) [kg/m^2]
    real(r8),intent(out) :: snow_depth(host%ncolumns_g)! snow depth [m]
    real(r8),intent(out) :: T_skin(host%ncolumns_g) ! skin temperature [K]
    real(r8),intent(out) :: T_veg(host%ncolumns_g) ! leaf temperature [K]
    real(r8),intent(out) :: T_soil(host%ncells_g)   ! soil temperature [K]
    call clm_to_host_diagnostics(host, clm, swe, snow_depth, canopy_storage, T_skin, T_veg, T_soil)
  end subroutine clm_to_ats_diagnostics

end module ats_clm
