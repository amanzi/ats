!Memory and interface mangling functions for the ATS-EcoSIM coupler.

module bgc_fortran_memory_mod


  use BGCContainers_module, only : BGCSizes,BGCProperties,&
           BGCState,BGCAuxiliaryData
  use iso_c_binding, only : c_funptr

  implicit none

  type, bind(C) :: BGCInterface

    type(c_funptr) :: DataTest
    type(c_funptr) :: Setup
    type(c_funptr) :: Shutdown
    type(c_funptr) :: Advance

  end type BGCInterface

  type, public :: BGCFortranInterface
    type(BGCInterface) :: c_interface

  contains
    procedure, public :: CreateInterface => Create_Fortran_BGC_Interface
    procedure, public :: DataTest => BGC_Fortran_DataTest
    procedure, public :: Setup  =>  BGC_Fortran_Setup
    procedure, public :: Shutdown  => BGC_Fortran_Shutdown
    procedure, public :: Advance => BGC_Fortran_Advance
  end type BGCFortranInterface

  interface
    subroutine CreateBGCInterface(engine_name, bgc_interface) bind(C, name='CreateBGCInterface')
      use iso_c_binding, only: c_char
      IMPORT
      implicit none
      character(kind=c_char) :: engine_name(*)
      type(BGCInterface)      :: bgc_interface
    end subroutine
  end interface

  ! Memory allocation subroutines

  interface
    subroutine AllocateBGCState(sizes, state, ncells_per_col_, num_components, num_columns) bind(C, name='AllocateBGCState')
      use BGCContainers_module, only : BGCSizes, BGCState
      use, intrinsic :: iso_c_binding, only: c_int
      implicit none
      type(BGCSizes) :: sizes
      type(BGCState) :: state
      integer(c_int),VALUE :: ncells_per_col_
      integer(c_int),VALUE :: num_components
      integer(c_int),VALUE :: num_columns
    end subroutine
  end interface
  interface
    subroutine FreeBGCState(state) bind(C, name='FreeBGCState')
      use BGCContainers_module, only : BGCState
      implicit none
      type(BGCState) :: state
    end subroutine
  end interface

  interface
    subroutine AllocateBGCProperties(sizes, properties, ncells_per_col_, num_columns) bind(C, name='AllocateBGCProperties')
      use BGCContainers_module, only : BGCSizes, BGCProperties
      use, intrinsic :: iso_c_binding, only: c_int
      implicit none
      type(BGCSizes) :: sizes
      type(BGCProperties) :: properties
      integer(c_int),VALUE :: ncells_per_col_
      integer(c_int),VALUE :: num_columns
    end subroutine
  end interface
  interface
    subroutine FreeBGCProperties(properties) bind(C, name='FreeBGCProperties')
      use BGCContainers_module, only : BGCProperties
      implicit none
      type(BGCProperties) :: properties
    end subroutine
  end interface

  ! The following subroutines are methods of the engine itself

  interface
    subroutine DataTest() bind(C)
      use, intrinsic :: iso_c_binding
      IMPORT
      implicit none
    end subroutine
  end interface


  interface
    subroutine Setup(properties, state, sizes, num_iterations, num_columns, ncells_per_col_) bind(C)

      use, intrinsic :: iso_c_binding, only: c_char, c_bool, c_ptr, c_int
      use BGCContainers_module, only : BGCSizes,BGCProperties,&
               BGCState
      IMPORT
      implicit none

      integer(c_int),VALUE :: num_iterations
      integer(c_int),VALUE :: num_columns
      integer(c_int),VALUE :: ncells_per_col_

      type(BGCProperties) :: properties
      type(BGCState) :: state
      type(BGCSizes) :: sizes

    end subroutine
  end interface

  !gracefully shutdown the engine, cleanup memory
  interface
    subroutine Shutdown() bind(C)
      use, intrinsic :: iso_c_binding, only : c_ptr

      implicit none

    end subroutine
  end interface

  ! take one (or more?) reaction steps in operator split mode
  interface
    subroutine Advance(delta_t, properties, state, sizes, num_iterations, num_columns) bind(C)
      use, intrinsic :: iso_c_binding, only : c_ptr, c_double, c_int
      use BGCContainers_module, only : BGCSizes,BGCProperties,&
               BGCState
      implicit none

      real(c_double),VALUE :: delta_t
      integer(c_int),VALUE :: num_iterations
      integer(c_int),VALUE :: num_columns

      type(BGCProperties) :: properties
      type(BGCState) :: state
      type(BGCSizes) :: sizes
    end subroutine
  end interface

  contains

  subroutine BGC_Fortran_DataTest(this)
    use, intrinsic :: iso_c_binding

    implicit none
    class(BGCFortranInterface) :: this
    procedure(DataTest), pointer :: engine_DataTest

    call c_f_procpointer(this%c_interface%Setup,engine_DataTest)
    call engine_DataTest()

  end subroutine BGC_Fortran_DataTest

  subroutine BGC_Fortran_Setup(this, properties, state, sizes, num_iterations,&
                               num_columns, ncells_per_col_)
    use, intrinsic :: iso_c_binding, only : c_ptr, c_int, c_double, c_f_procpointer
    use BGCContainers_module, only : BGCSizes, BGCProperties,&
             BGCState

    implicit none
    class(BGCFortranInterface) :: this

    real(c_double) :: delta_t
    integer(c_int) :: num_columns
    integer(c_int) :: num_iterations
    integer(c_int) :: ncells_per_col_
    type(BGCProperties) :: properties
    type(BGCState) :: state
    type(BGCSizes) :: sizes

    procedure(Setup), pointer :: engine_Setup

    call c_f_procpointer(this%c_interface%Setup,engine_Setup)
    call engine_Setup(properties, state, sizes, num_iterations, &
                      num_columns, ncells_per_col_)

  end subroutine BGC_Fortran_Setup

  subroutine BGC_Fortran_Shutdown(this)
    use, intrinsic :: iso_c_binding, only : c_ptr,c_f_procpointer

    implicit none
    class(BGCFortranInterface) :: this
    procedure(Shutdown), pointer :: engine_Shutdown

    call c_f_procpointer(this%c_interface%Shutdown,engine_Shutdown)
    call engine_Shutdown()

  end subroutine BGC_Fortran_Shutdown

  subroutine BGC_Fortran_Advance(this, delta_t, properties, state, sizes, num_iterations, num_columns)
    use, intrinsic :: iso_c_binding, only : c_ptr, c_int, c_double, c_f_procpointer
    use BGCContainers_module, only : BGCSizes, BGCProperties,&
             BGCState

    implicit none
    class(BGCFortranInterface) :: this

    real(c_double) :: delta_t
    integer(c_int) :: num_columns
    integer(c_int) :: num_iterations
    type(BGCProperties) :: properties
    type(BGCState) :: state
    type(BGCSizes) :: sizes

    procedure(Advance), pointer :: engine_Advance

    call c_f_procpointer(this%c_interface%Advance,engine_Advance)
    call engine_Advance(delta_t, properties, state, sizes, num_iterations, num_columns)
  end subroutine BGC_Fortran_Advance

  subroutine Create_Fortran_BGC_Interface(this,engine_name)
    use BGCContainers_module, only :kBGCMaxStringLength
    use iso_c_binding, only: c_char,c_null_char
    implicit none
    class(BGCFortranInterface) :: this
    character(kind=c_char,len=kBGCMaxStringLength) :: engine_name

    call CreateBGCInterface(trim(engine_name)//C_NULL_CHAR, this%c_interface)

  end subroutine

end module bgc_fortran_memory_mod
