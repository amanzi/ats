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

set(elm_ats_src_files
  ${ATS_SOURCE_DIR}/src/executables/coordinator.cc
  ${ATS_SOURCE_DIR}/src/executables/ats_mesh_factory.cc
  elm_ats_driver.cc
  elm_ats_api.cc
  )

set(elm_ats_inc_files
  ${ATS_SOURCE_DIR}/src/executables/coordinator.hh
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
  set_target_properties(elm_ats PROPERTIES LINK_FLAGS " -ld_classic -Wl,-undefined,dynamic_lookup")
endif()

# ## build static lib
# add_amanzi_library(elm_ats_static
#   STATIC
#   SOURCE ${elm_ats_src_files}
#   HEADERS ${elm_ats_inc_files}
#   LINK_LIBS ${fates_link_libs} ${tpl_link_libs} ${ats_link_libs} ${amanzi_link_libs})

# re-enable later
##add_subdirectory(test)
