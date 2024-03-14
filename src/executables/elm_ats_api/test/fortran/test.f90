
program elm_test
    use libats
    implicit none
    include 'mpif.h'
    type(ats) :: ats_driver
    Character(len = 200) :: infile_name
    integer :: ierror, i

    ! dummy data
    integer, parameter :: ncol = 5
    integer, parameter :: ncell = 15
    double precision, dimension(ncol*ncell) :: soil_pres
    double precision, dimension(ncol*ncell) :: satur
    double precision, dimension(ncol) :: tran
    double precision, dimension(ncol) :: infil
    double precision, dimension(ncol) :: evap
    integer :: ncols_local, ncols_global, ncells_per_col
    double precision :: time
    character (len = 40) :: test_prefix

    infil(:) = 10.0
    evap(:) = 3.0
    tran(:) = 6.0

    do i=1,ncol
      tran(i) = (1.0 - (1.0/ncell)*i)
    end do

    call get_command_argument(1, infile_name)

    ! remove prefix if sent from test suite
    if (infile_name(1:11) == '--xml_file=') then
      infile_name = trim(adjustl(infile_name(12:)))
    end if

    time = 0.0

    call MPI_INIT(ierror)

    ! create an ATS driver object
    ats_driver = ats(MPI_COMM_WORLD, infile_name)

    ! call ATS methods
    call ats_driver%setup()
    call ats_driver%initialize(time, satur, soil_pres)
    !call ats_driver%set_sources(infil, evap, tran)
    call ats_driver%advance_test()

    ! don't need to call ats_driver%delete now that final is used

    print*,"DONE WITH ELM-ATS FORTRAN TEST"

    call MPI_FINALIZE(ierror)

end program elm_test

