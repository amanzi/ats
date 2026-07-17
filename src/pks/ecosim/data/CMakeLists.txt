# -*- mode: cmake -*-
#
#  ATS
#    Data management programs for EcoSIM_ATS
#

# get_property(AMANZI_TPLS_DIR GLOBAL PROPERTY AMANZI_TPLS_DIR)

set(ECOSIM_INSTALL_PREFIX ${ECOSIM_DIR})
set(ECOSIM_INCLUDE_DIRS ${ECOSIM_DIR}/include)
set(ECOSIM_LIB_DIR ${ECOSIM_DIR}/lib)

message("In ATS-EcoSIM data CMakeLists:")
message("ECOSIM_DIR:" ${ECOSIM_DIR})
message("ECOSIM_INSTALL_PREFIX:" ${ECOSIM_INSTALL_PREFIX})
message("ECOSIM_INCLUDE_DIRS:" ${ECOSIM_INCLUDE_DIRS})
message("ECOSIM_LIB_DIR:" ${ECOSIM_LIB_DIR})

message("NetCDF_INCLUDE_DIRS: ${NetCDF_INCLUDE_DIRS}")
message("NetCDF_C_LIBRARIES: ${NetCDF_C_LIBRARIES}")

include_directories(${ECOSIM_INCLUDE_DIRS})

set(ats_ecosim_data_src_files
    BGC_memory.cc
    BGC_containers.cc
    bgc_fortran_memory_mod.F90
)

set(ats_ecosim_data_inc_files
  BGC_containers.hh
  BGC_memory.hh
)

file(GLOB ECOSIM_LIBRARIES
 ${ECOSIM_LIB_DIR}/*
)

set(ats_ecosim_data_link_libs
  ${Teuchos_LIBRARIES}
  ${Epetra_LIBRARIES}
  ${ECOSIM_LIBRARIES}
  error_handling
  atk
  mesh
  data_structures
  whetstone
  operators
  solvers
  time_integration
  state
  pks
  chemistry_pk
  ats_pks
  ats_eos
  ats_operators
  ${NetCDF_C_LIBRARIES}
)

message(STATUS "ats_ecosim_data_link_libs: ${ats_ecosim_data_link_libs}")

# make the library
add_amanzi_library(ats_ecosim_data
                   SOURCE ${ats_ecosim_data_src_files}
                   HEADERS ${ats_ecosim_data_inc_files}
                   LINK_LIBS ${ats_ecosim_data_link_libs}
                   )

generate_evaluators_registration_header(
    HEADERFILE ats_ecosim_data_registration.hh
    LISTNAME   ATS_ECOSIM_DATA_REG
    INSTALL    True
    )
