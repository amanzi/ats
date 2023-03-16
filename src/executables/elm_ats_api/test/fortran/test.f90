
program elm_test
    use libats
    implicit none
    include 'mpif.h'
    type(ats) :: ats_driver
    Character(len = 200) :: infile_name
    integer :: ierror, i

    ! dummy data
    ! 1 column, 100 cells
    integer, parameter :: ncol = 1
    integer, parameter :: ncell = 100
    double precision, dimension(ncol) :: infil
    double precision, dimension(ncol) :: evap
    double precision, dimension(ncol) :: surf_pres
    double precision, dimension(ncol) :: elev
    double precision, dimension(ncol) :: surf_area_m2
    double precision, dimension(ncol) :: lat
    double precision, dimension(ncol) :: lon

    double precision, dimension(ncell) :: dz
    double precision, dimension(ncell) :: depth
    double precision, dimension(ncell) :: soil_pres
    double precision, dimension(ncell) :: satur
    double precision, dimension(ncell) :: tran

    integer :: ncols_local, ncols_global, ncells_per_col
    character (len = 40) :: test_prefix

    infil(:) = 10.0
    evap(:) = 3.0

    do i=1,ncell
      tran(i) = (1.0 - (1.0/ncell)*i)
    end do

    call get_command_argument(1, infile_name)

    ! remove prefix if sent from test suite
    if (infile_name(1:11) == '--xml_file=') then
      infile_name = trim(adjustl(infile_name(12:)))
    end if

    call MPI_INIT(ierror)

    ! create an ATS driver object
    ats_driver = ats(MPI_COMM_WORLD, infile_name)

    ! call ATS methods
    call ats_driver%setup()
    !!call ats_driver%get_mesh_info(ncols_local, ncols_global, ncells_per_col, dz, depth, elev, surf_area_m2, lat, lon)
    call ats_driver%initialize()
    !!call ats_driver%set_sources(infil, evap, tran, ncol, ncell)
    call ats_driver%advance_test()
    !!call ats_driver%get_waterstate(surf_pres, soil_pres, satur, ncol, ncell)

    ! don't need to call ats_driver%delete now that final is used

    print*,"DONE WITH ELM-ATS FORTRAN TEST"

    call MPI_FINALIZE(ierror)

end program elm_test

