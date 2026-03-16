# -*- mode: cmake -*-
#
#  ATS
#    Data management programs for EcoSIM_ATS
#

get_property(AMANZI_TPLS_DIR GLOBAL PROPERTY AMANZI_TPLS_DIR)

set(ECOSIM_INSTALL_PREFIX ${ECOSIM_DIR}/ecosim)
set(ECOSIM_LIB_LOCATION ${ECOSIM_DIR}/ecosim/local/lib)
set(ECOSIM_BUILD_PREFIX ${ECOSIM_DIR}/ecosim/build)
set(NETCDF_LIB ${ECOSIM_DIR}/lib)

message("In ATS-EcoSIM data CMakeLists:")
message("ECOSIM_DIR:" ${ECOSIM_DIR})
message("ECOSIM_INSTALL_PREFIX:" ${ECOSIM_INSTALL_PREFIX})
message("ECOSIM_LIB_LOCATION:" ${ECOSIM_LIB_LOCATION})
message("ECOSIM_BUILD_PREFIX:" ${ECOSIM_BUILD_PREFIX})

message("At include_directories in data")
include_directories(${ECOSIM_BUILD_PREFIX}/f90src/Utils/)
include_directories(${ECOSIM_BUILD_PREFIX}/f90src/Minimath/)
include_directories(${ECOSIM_BUILD_PREFIX}/f90src/Modelconfig/)
include_directories(${ECOSIM_BUILD_PREFIX}/f90src/Modelforc/)
include_directories(${ECOSIM_BUILD_PREFIX}/f90src/Mesh/)
include_directories(${ECOSIM_BUILD_PREFIX}/f90src/Modelpars/)
include_directories(${ECOSIM_BUILD_PREFIX}/f90src/Balances/)
include_directories(${ECOSIM_BUILD_PREFIX}/f90src/Ecosim_datatype/)
include_directories(${ECOSIM_BUILD_PREFIX}/f90src/HydroTherm/SoilPhys/)
include_directories(${ECOSIM_BUILD_PREFIX}/f90src/HydroTherm/SurfPhys/)
include_directories(${ECOSIM_BUILD_PREFIX}/f90src/HydroTherm/SnowPhys/)
include_directories(${ECOSIM_BUILD_PREFIX}/f90src/HydroTherm/PhysData/)
include_directories(${ECOSIM_BUILD_PREFIX}/f90src/Ecosim_mods/)
include_directories(${ECOSIM_BUILD_PREFIX}/f90src/Prescribed_pheno/)
include_directories(${ECOSIM_BUILD_PREFIX}/f90src/Plant_bgc/)
include_directories(${ECOSIM_BUILD_PREFIX}/f90src/APIs/)
include_directories(${ECOSIM_BUILD_PREFIX}/f90src/APIData/)

set(ats_ecosim_data_src_files
    BGC_memory.cc
    BGC_containers.cc
    bgc_fortran_memory_mod.F90
    ${ECOSIM_INSTALL_PREFIX}/f90src/ATSUtils/BGC_containers.F90
    ${ECOSIM_INSTALL_PREFIX}/f90src/ATSUtils/BGC_containers.F90
    ${ECOSIM_INSTALL_PREFIX}/f90src/ATSUtils/ATSCPLMod.F90
    ${ECOSIM_INSTALL_PREFIX}/f90src/ATSUtils/ATSEcoSIMInitMod.F90
    ${ECOSIM_INSTALL_PREFIX}/f90src/ATSUtils/SharedDataMod.F90
    ${ECOSIM_INSTALL_PREFIX}/f90src/ATSUtils/c_f_interface_module.F90
    ${ECOSIM_INSTALL_PREFIX}/f90src/ATSUtils/ATSEcoSIMAdvanceMod.F90
    ${ECOSIM_INSTALL_PREFIX}/f90src/ATSUtils/ATSUtilsMod.F90
)

set(ats_ecosim_data_inc_files
  BGC_containers.hh
  BGC_memory.hh
)

file(GLOB ECOSIM_LIBRARIES
 ${ECOSIM_LIB_LOCATION}/*.a
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
  ${NETCDF_LIB}/libnetcdf.so
  ${NETCDF_LIB}/libnetcdf.so.22
  ${NETCDF_LIB}/libnetcdff.so
  ${NETCDF_LIB}/libnetcdff.so.7
  ${NETCDF_LIB}/libnetcdff.so.7.1.0
  gfortran
)

message(STATUS "ats_ecosim_data_link_libs: ${ats_ecosim_data_link_libs}")

# make the library
add_amanzi_library(ats_ecosim_data
                   SOURCE ${ats_ecosim_data_src_files}
                   HEADERS ${ats_ecosim_data_inc_files}
     LINK_LIBS ${ats_ecosim_data_link_libs})

generate_evaluators_registration_header(
    HEADERFILE ats_ecosim_data_registration.hh
    LISTNAME   ATS_ECOSIM_DATA_REG
    INSTALL    True
    )
