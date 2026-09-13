module elm_water_balance_kernel

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Dependency-free kernel for the ELM column water-balance error (errh2o).
  !
  ! This module is deliberately self-contained (it uses only iso_c_binding)
  ! so that the *same source file* can be compiled into both:
  !   (1) the ELM/E3SM build, where ColWaterBalanceCheck calls
  !       elm_water_balance_error() so ELM's balance check has a single
  !       source of truth, and
  !   (2) the ATS elm_ats_api library, where the C-bound entry point
  !       elm_ats_water_balance_error_c() is called across the coupling
  !       boundary to gate ATS step acceptance on the ELM balance.
  !
  ! Keeping one file compiled into both codes is what lets ATS reuse ELM's
  ! balance formula without a second, drift-prone re-implementation.
  !
  ! Sign / unit conventions (identical to BalanceCheckMod.F90):
  !   - all storages in kg m-2 (mm H2O), all fluxes in mm H2O s-1,
  !     dtime in s; result errh2o in mm H2O.
  !   - sources (precip, irrigation, flood, from_uphill, floodplain drain
  !     in) enter with +, evaporation/transpiration and all runoff/drainage
  !     terms are sinks (subtracted).
  !-----------------------------------------------------------------------
  use iso_c_binding, only : c_double

  implicit none
  private

  ! kind used throughout; c_double matches ELM's shr_kind_r8 (8-byte real)
  integer, parameter :: wb_r8 = c_double

  public :: elm_water_balance_error       ! pure Fortran formula (used by ELM)
  public :: elm_ats_water_balance_error_c ! C-bound entry (used by ATS)

contains

  !-----------------------------------------------------------------------
  ! Full ELM column water-balance error, one column.
  ! This is the arithmetic core formerly inlined in ColWaterBalanceCheck
  ! (BalanceCheckMod.F90).  errh2o = dStorage - net_flux * dtime.
  !-----------------------------------------------------------------------
  pure function elm_water_balance_error( &
       endwb, begwb, forc_rain, forc_snow, qflx_floodc, qflx_from_uphill, &
       qflx_surf_irrig, qflx_over_supply, qflx_evap_tot, qflx_surf, &
       qflx_h2osfc_surf, qflx_to_downhill, qflx_qrgwl, qflx_drain, &
       qflx_drain_perched, qflx_snwcp_ice, qflx_ice_runoff_xs, qflx_lateral, &
       qflx_h2orof_drain, qflx_lnd2ocn, qflx_h2oocn_drain, dtime) result(errh2o)
    real(wb_r8), intent(in) :: endwb, begwb
    real(wb_r8), intent(in) :: forc_rain, forc_snow, qflx_floodc, qflx_from_uphill
    real(wb_r8), intent(in) :: qflx_surf_irrig, qflx_over_supply, qflx_evap_tot, qflx_surf
    real(wb_r8), intent(in) :: qflx_h2osfc_surf, qflx_to_downhill, qflx_qrgwl, qflx_drain
    real(wb_r8), intent(in) :: qflx_drain_perched, qflx_snwcp_ice, qflx_ice_runoff_xs, qflx_lateral
    real(wb_r8), intent(in) :: qflx_h2orof_drain, qflx_lnd2ocn, qflx_h2oocn_drain, dtime
    real(wb_r8) :: errh2o

    errh2o = endwb - begwb &
         - (forc_rain + forc_snow + qflx_floodc + qflx_from_uphill &
         + qflx_surf_irrig + qflx_over_supply &
         - qflx_evap_tot - qflx_surf - qflx_h2osfc_surf - qflx_to_downhill &
         - qflx_qrgwl - qflx_drain - qflx_drain_perched - qflx_snwcp_ice - qflx_ice_runoff_xs &
         - qflx_lateral + qflx_h2orof_drain - qflx_lnd2ocn + qflx_h2oocn_drain) * dtime
  end function elm_water_balance_error

  !-----------------------------------------------------------------------
  ! C-bound entry for ATS.  ATS solves only surface + subsurface water, so
  ! it supplies the subset of ELM flux terms it controls and the routine
  ! passes zero for the ELM-only terms (snow capping, glacier, irrigation,
  ! flood, lateral hillslope, ocean exchange).
  !
  ! Mapping (ATS term -> ELM term in the full formula above):
  !   source    -> forc_rain          (net water into soil top)
  !   evap+tran -> qflx_evap_tot
  !   runoff    -> qflx_surf           (surface runoff)
  !   baseflow  -> qflx_drain          (subsurface runoff)
  !   dStorage  -> endwb - begwb
  ! All arguments are passed by reference (Fortran convention); the C++
  ! caller passes double* .  Result returned by value as a C double.
  !-----------------------------------------------------------------------
  function elm_ats_water_balance_error_c(endwb, begwb, source, evap, tran, &
       baseflow, runoff, dtime) result(errh2o) &
       bind(C, name="elm_ats_water_balance_error_c")
    real(c_double), intent(in) :: endwb, begwb
    real(c_double), intent(in) :: source, evap, tran, baseflow, runoff, dtime
    real(c_double) :: errh2o
    real(c_double), parameter :: z = 0.0_c_double

    errh2o = elm_water_balance_error( &
         endwb, begwb, &
         source, z, z, z, &        ! forc_rain=source, forc_snow, qflx_floodc, qflx_from_uphill
         z, z, evap + tran, runoff, & ! qflx_surf_irrig, qflx_over_supply, qflx_evap_tot, qflx_surf
         z, z, z, baseflow, &      ! qflx_h2osfc_surf, qflx_to_downhill, qflx_qrgwl, qflx_drain
         z, z, z, z, &             ! qflx_drain_perched, qflx_snwcp_ice, qflx_ice_runoff_xs, qflx_lateral
         z, z, z, dtime)           ! qflx_h2orof_drain, qflx_lnd2ocn, qflx_h2oocn_drain, dtime
  end function elm_ats_water_balance_error_c

end module elm_water_balance_kernel
