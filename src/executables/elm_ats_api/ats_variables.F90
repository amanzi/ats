module ats_variables
  use, intrinsic :: iso_c_binding
  implicit none

  type :: ats_var_id_type
     integer(c_int), parameter :: NULL = 0
     integer(c_int), parameter :: BASE_POROSITY = 1
     integer(c_int), parameter :: HYDRAULIC_CONDUCTIVITY = 2
     integer(c_int), parameter :: CLAPP_HORNBERGER_B = 3
     integer(c_int), parameter :: CLAPP_HORNBERGER_PSI_SAT = 4
     integer(c_int), parameter :: RESIDUAL_SATURATION = 5
     integer(c_int), parameter :: EFFECTIVE_POROSITY = 6
     integer(c_int), parameter :: ROOT_FRACTION = 7
     integer(c_int), parameter :: SURFACE_WATER_CONTENT = 8
     integer(c_int), parameter :: WATER_CONTENT = 9
     integer(c_int), parameter :: PRESSURE = 10
     integer(c_int), parameter :: GROSS_SURFACE_WATER_SOURCE = 11
     integer(c_int), parameter :: POTENTIAL_EVAPORATION = 12
     integer(c_int), parameter :: POTENTIAL_TRANSPIRATION = 13
     integer(c_int), parameter :: EVAPORATION = 14
     integer(c_int), parameter :: COLUMN_TRANSPIRATION = 15
     integer(c_int), parameter :: RUNOFF = 16
     integer(c_int), parameter :: BASEFLOW = 17
     integer(c_int), parameter :: ELEVATION = 18
     integer(c_int), parameter :: TIME = 19
  end type ats_var_id_type

  type(ats_var_id_type), parameter :: ats_var_id = ats_var_id_type()
  
end module ats_variables

