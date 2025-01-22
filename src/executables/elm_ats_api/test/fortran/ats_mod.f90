
!! Defines a Fortran type that can store a pointer to ATS
!! Provides implementations of methods declared in ats.f90
module libats
    use iso_c_binding
    private
    public :: ats
    include "ats.f90"

    type ats
        private
        type(c_ptr) :: ptr ! pointer to ats class
    contains
        final :: ats_delete 
        procedure :: setup => ats_setup
        procedure :: initialize => ats_initialize
        procedure :: advance => ats_advance
        procedure :: advance_test => ats_advance_test
        procedure :: set_sources => ats_set_sources
        procedure :: get_waterstate => ats_get_waterstate
        procedure :: get_mesh_info => ats_get_mesh_info
    end type

    ! constructor
    interface ats
        procedure ats_create
    end interface

!! Definitions of Fortran interface methods
!! These methods act like wrappers around their C bindings.
contains
    
    !! Allocate and construct a new ATS object
    !! Returns a pointer to new ATS object
    function ats_create(comm, infile)
        implicit none
        type(ats) :: ats_create
        integer, intent(in) :: comm
        character(len=*), intent(in) :: infile
        character(len=1, kind=C_CHAR) :: c_str_infile(len_trim(infile) + 1)
        integer :: n_char, i

        !! basic parsing
        n_char = len_trim(infile)
        do i = 1, n_char
            c_str_infile(i) = infile(i:i)
        end do
        c_str_infile(n_char + 1) = C_NULL_CHAR

        ats_create%ptr = ats_create_c(comm, c_str_infile)
    end function

    !! Destroy ATS object
    subroutine ats_delete(this)
        implicit none
        type(ats) :: this
        call ats_delete_c(this%ptr)
    end subroutine

    !! Setup ATS for LSM simulation
    subroutine ats_setup(this)
        implicit none
        class(ats) :: this
       call ats_setup_c(this%ptr)
    end subroutine

    !! Initialize ATS for this simulation
    subroutine ats_initialize(this)
        implicit none
        class(ats) :: this
        call ats_initialize_c(this%ptr)
    end subroutine

    !! Advance ATS by dt
    subroutine ats_advance(this, dt)
        implicit none
        class(ats) :: this
        double precision, intent(in) :: dt
        call ats_advance_c(this%ptr, dt)
    end subroutine

    !! Temporary for testing - implements time loop and calls ats_advance()
    subroutine ats_advance_test(this)
        implicit none
        class(ats) :: this
        call ats_advance_test_c(this%ptr)
    end subroutine

    !! WIP
    !! Sets ATS coupling variables infiltration, evaporation, and transpiration
    subroutine ats_set_sources(this, soil_infil, soil_evap, root_trans, ncols, ncells)
        implicit none
        class(ats) :: this
        double precision, dimension(*), intent(in) :: soil_infil
        double precision, dimension(*), intent(in) :: soil_evap
        double precision, dimension(*), intent(in) :: root_trans
        integer, intent(in) :: ncols
        integer, intent(in) :: ncells
        call ats_set_sources_c(this%ptr, soil_infil, soil_evap, root_trans, ncols, ncells)
    end subroutine

    !! WIP
    subroutine ats_get_waterstate(this, surf_pres, soil_pres, satur, ncols, ncells)
        implicit none
        class(ats) :: this
        double precision, dimension(*), intent(in) :: surf_pres
        double precision, dimension(*), intent(in) :: soil_pres
        double precision, dimension(*), intent(in) :: satur
        integer, intent(in) :: ncols
        integer, intent(in) :: ncells
        call ats_get_waterstate_c(this%ptr, surf_pres, soil_pres, satur, ncols, ncells)
    end subroutine

    !! WIP
    subroutine ats_get_mesh_info(this, ncols_local, ncols_global, ncells_per_col, dz, depth, elev, surf_area_m2, lat, lon)
        implicit none
        class(ats) :: this
        double precision, dimension(*), intent(in) :: dz
        double precision, dimension(*), intent(in) :: depth
        double precision, dimension(*), intent(in) :: elev
        double precision, dimension(*), intent(in) :: surf_area_m2
        double precision, dimension(*), intent(in) :: lat
        double precision, dimension(*), intent(in) :: lon
        integer, intent(in) :: ncols_local
        integer, intent(in) :: ncols_global
        integer, intent(in) :: ncells_per_col
        call ats_get_mesh_info_c(this%ptr, ncols_local, ncols_global, ncells_per_col, dz, depth, elev, surf_area_m2, lat, lon)
    end subroutine

end module
