
program elm_test
    use libats
    implicit none
    include 'mpif.h'
    type(ats) :: ats_driver
    Character(len = 200) :: infile_name
    integer :: ierror, i

    ! dummy data
    ! 1 column, 100 cells
    double precision, dimension(1) :: infil
    double precision, dimension(1) :: evap
    double precision, dimension(100) :: tran

    double precision, dimension(1) :: surf_pres
    double precision, dimension(100) :: soil_pres
    double precision, dimension(100) :: satur

    infil(1) = 10.0
    evap(1) = 3.0

    do i=1,100
      tran(i) = (1.0 - 0.01*i)
    end do

    call get_command_argument(1, infile_name)

    call MPI_INIT(ierror)

    ! create an ATS driver object
    ats_driver = ats()

    ! call ATS methods
    call ats_driver%setup(MPI_COMM_WORLD, infile_name)
    call ats_driver%initialize()
    call ats_driver%set_sources(infil, evap, tran, 1, 100)
    call ats_driver%advance_test()

    call ats_driver%get_waterstate(surf_pres, soil_pres, satur, 1, 100)

    ! don't need to call anymore now that final is used
    !call ats_driver%delete

    print*,"DONE WITH ELM-ATS FORTRAN TEST"

    call MPI_FINALIZE(ierror)

end program elm_test

