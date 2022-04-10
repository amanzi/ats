
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
    end type

    ! wow fortran constructor wow
    interface ats
        procedure ats_create
    end interface

contains
    
    function ats_create()
        implicit none
        type(ats) :: ats_create
        ats_create%ptr = ats_create_c()
    end function

    subroutine ats_delete(this)
        implicit none
        type(ats) :: this
        call ats_delete_c(this%ptr)
    end subroutine

    integer function ats_setup(this, infile)
        implicit none
        class(ats) :: this
        character(len=*), intent(in) :: infile
        character(len=1, kind=C_CHAR) :: c_str_infile(len_trim(infile) + 1)
        integer :: n_char, i

        n_char = len_trim(infile)
        do i = 1, n_char
            c_str_infile(i) = infile(i:i)
        end do
        c_str_infile(n_char + 1) = C_NULL_CHAR

        ats_setup = ats_setup_c(this%ptr, c_str_infile)
    end function

    subroutine ats_initialize(this)
        implicit none
        class(ats) :: this
        call ats_initialize_c(this%ptr)
    end subroutine

    subroutine ats_advance(this, dt)
        implicit none
        class(ats) :: this
        double precision, intent(in) :: dt
        call ats_advance_c(this%ptr, dt)
    end subroutine

    subroutine ats_advance_test(this)
        implicit none
        class(ats) :: this
        call ats_advance_test_c(this%ptr)
    end subroutine

end module
