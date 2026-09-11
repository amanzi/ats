# -*- mode: cmake -*-
#
#  ELM_ATS_API
#    library
#
# builds both dynamic and static libs for now
# need to add build logic to make configurable

# verify the compatibility of the C/Fortran and C++/Fortran compilers
# Amanzi probably does this somewhere upstream, should check
include(FortranCInterface)
FortranCInterface_VERIFY(CXX)

include_directories(${ATS_SOURCE_DIR}/src/executables)

# Shared ELM water-balance kernel: the *same* dependency-free Fortran source
# that ELM compiles into BalanceCheckMod is compiled here so ATS calls ELM's
# balance formula (via bind(C)) rather than duplicating it.  Resolve its path
# from E3SM_SRC_DIR (set per the project environment; see CLAUDE.md).
if (DEFINED ENV{E3SM_SRC_DIR})
  set(_e3sm_src_dir "$ENV{E3SM_SRC_DIR}")
else()
  set(_e3sm_src_dir "${ATS_SOURCE_DIR}/../../../../../E3SM")
endif()
set(elm_water_balance_kernel_src
  "${_e3sm_src_dir}/components/elm/src/biogeophys/elm_water_balance_kernel.F90")
if (NOT EXISTS "${elm_water_balance_kernel_src}")
  message(FATAL_ERROR
    "elm_ats_api: shared ELM water-balance kernel not found at "
    "${elm_water_balance_kernel_src}. Set E3SM_SRC_DIR to the E3SM checkout.")
endif()

set(elm_ats_src_files
  ${ATS_SOURCE_DIR}/src/executables/driver.cc
  ${ATS_SOURCE_DIR}/src/executables/ats_mesh_factory.cc
  ${elm_water_balance_kernel_src}
  elm_ats_driver.cc
  elm_ats_api.cc
  )

set(elm_ats_inc_files
  ${ATS_SOURCE_DIR}/src/executables/driver.hh
  ${ATS_SOURCE_DIR}/src/executables/ats_mesh_factory.hh
  elm_ats_driver.hh
  elm_ats_api.h
  )

## build shared lib
add_amanzi_library(elm_ats
  SOURCE ${elm_ats_src_files}
  HEADERS ${elm_ats_inc_files}
  LINK_LIBS ${fates_link_libs} ${tpl_link_libs} ${ats_link_libs} ${amanzi_link_libs})
if (APPLE AND BUILD_SHARED_LIBS)
  set_target_properties(elm_ats PROPERTIES LINK_FLAGS "-Wl,-undefined,dynamic_lookup")
endif()

# ## build static lib
# add_amanzi_library(elm_ats_static
#   STATIC
#   SOURCE ${elm_ats_src_files}
#   HEADERS ${elm_ats_inc_files}
#   LINK_LIBS ${fates_link_libs} ${tpl_link_libs} ${ats_link_libs} ${amanzi_link_libs})

# re-enable later
##add_subdirectory(test)
