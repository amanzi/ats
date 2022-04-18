
program elm_test
    use libats
    implicit none
    include 'mpif.h'
    type(ats) :: ats_driver
    Character(len = 200) :: infile_name
    integer :: ierror

    call get_command_argument(1, infile_name)

    call MPI_INIT(ierror)

    ! create an ATS driver object
    ats_driver = ats()

    ! call ATS methods
    call ats_driver%setup(MPI_COMM_WORLD, infile_name)
    call ats_driver%initialize()
    call ats_driver%advance_test()

    ! don't need to call anymore now that final is used
    !call ats_driver%delete

    call MPI_FINALIZE(ierror)

end program elm_test

