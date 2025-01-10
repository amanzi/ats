
!! Defines a Fortran interface to ATS library
!! Provides Fortran method declarations and bindings to ATS C API functions
!! The implementation of the methods declared here can be found in ats_mod.f90
interface

    function ats_create_c(comm, infile) bind(c, name="ats_create")
        use iso_c_binding
        implicit none
        type(c_ptr) :: ats_create_c
        integer, intent(in) :: comm
        character(len=1, kind=C_CHAR), intent(in) :: infile(*)
    end function
    
    subroutine ats_delete_c(ats) bind(c, name="ats_delete")
        use iso_c_binding
        implicit none
        type(c_ptr) :: ats
    end subroutine
    
    subroutine ats_setup_c(ats) bind(c, name="ats_setup")
        use iso_c_binding
        implicit none
        type(c_ptr), value :: ats
    end subroutine
    
    subroutine ats_initialize_c(ats, time, satur, soil_pres) bind(c, name="ats_initialize")
        use iso_c_binding
        implicit none
        type(c_ptr), value :: ats
        real(c_double), intent(in) :: time
        real(c_double), dimension(*), intent(in) :: satur
        real(c_double), dimension(*), intent(in) :: soil_pres
    end subroutine
    
    subroutine ats_advance_test_c(ats) bind(c, name="ats_advance_test")
        use iso_c_binding
        implicit none
        type(c_ptr), value :: ats
    end subroutine
    
    subroutine ats_advance_c(ats, dt) bind(c, name="ats_advance")
        use iso_c_binding
        implicit none
        type(c_ptr), value :: ats
        real(c_double), intent(in) :: dt
    end subroutine

    subroutine ats_set_sources_c(ats, soil_infil, soil_evap, root_trans) &
        bind(c, name="ats_set_sources")
        use iso_c_binding
        implicit none
        type(c_ptr), value :: ats
        real(c_double), dimension(*), intent(in) :: soil_infil
        real(c_double), dimension(*), intent(in) :: soil_evap
        real(c_double), dimension(*), intent(in) :: root_trans
    end subroutine

    subroutine ats_get_waterstate_c(ats, surf_pres, soil_pres, satur, ncols, ncells) &
        bind(c, name="ats_get_waterstate")
        use iso_c_binding
        implicit none
        type(c_ptr), value :: ats
        real(c_double), dimension(*), intent(in) :: surf_pres
        real(c_double), dimension(*), intent(in) :: soil_pres
        real(c_double), dimension(*), intent(in) :: satur
        integer(c_int), intent(in) :: ncols
        integer(c_int), intent(in) :: ncells
    end subroutine

    subroutine ats_get_mesh_info_c(ats, ncols_local, ncols_global, ncells_per_col, dz, &
        depth, elev, surf_area_m2, lat, lon) bind(c, name="ats_get_mesh_info")
        use iso_c_binding
        implicit none
        type(c_ptr), value :: ats
        real(c_double), dimension(*), intent(in) :: dz
        real(c_double), dimension(*), intent(in) :: depth
        real(c_double), dimension(*), intent(in) :: elev
        real(c_double), dimension(*), intent(in) :: surf_area_m2
        real(c_double), dimension(*), intent(in) :: lat
        real(c_double), dimension(*), intent(in) :: lon
        integer(c_int), intent(in) :: ncols_local
        integer(c_int), intent(in) :: ncols_global
        integer(c_int), intent(in) :: ncells_per_col
    end subroutine

end interface
