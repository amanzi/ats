
!! test interface
interface

    function ats_create_c() bind(c, name="ats_create")
        use iso_c_binding
        implicit none
        type(c_ptr) :: ats_create_c
    end function
    
    subroutine ats_delete_c(ats) bind(c, name="ats_delete")
        use iso_c_binding
        implicit none
        type(c_ptr) :: ats
    end subroutine
    
    subroutine ats_setup_c(ats, comm, infile) bind(c, name="ats_setup")
        use iso_c_binding
        implicit none
        type(c_ptr), value :: ats
        integer, intent(in) :: comm
        character(len=1, kind=C_CHAR), intent(in) :: infile(*)
    end subroutine
    
    subroutine ats_initialize_c(ats) bind(c, name="ats_initialize")
        use iso_c_binding
        implicit none
        type(c_ptr), value :: ats
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

    subroutine ats_set_sources_c(ats, soil_infil, soil_evap, root_trans, ncols, ncells) bind(c, name="ats_set_sources")
        use iso_c_binding
        implicit none
        type(c_ptr), value :: ats
        real(c_double), dimension(*), intent(in) :: soil_infil
        real(c_double), dimension(*), intent(in) :: soil_evap
        real(c_double), dimension(*), intent(in) :: root_trans
        integer(c_int), intent(in) :: ncols
        integer(c_int), intent(in) :: ncells
    end subroutine

end interface
