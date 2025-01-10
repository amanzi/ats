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
include_directories(${ATS_SOURCE_DIR}/src/executables/elm_ats_api)

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
  set_target_properties(elm_ats PROPERTIES LINK_FLAGS "-Wl,-undefined,dynamic_lookup")
endif()


#
# These used to work, but need an inputfile that I'm not sure where it is... --ETC
# 
# if (BUILD_TESTS)

#   # Add UnitTest includes
#   include_directories(${UnitTest_INCLUDE_DIRS})
#   include_directories(${EXECUTABLE_SOURCE_DIR})

#   # test for ELM API
#   add_amanzi_test(elm_ats_cc_test elm_ats_cc_test
#     KIND int
#     SOURCE test/elm_ats_test.cc
#     LINK_LIBS elm_ats ${ats_link_libs} ${UnitTest_LIBRARIES})

#   # test for ELM API in fortran
#   add_amanzi_test(elm_ats_fort_test elm_ats_fort_test
#     KIND int
#     SOURCE test/fortran/ats_mod.f90 test/fortran/test.f90
#     LINK_LIBS elm_ats ${ats_link_libs} ${UnitTest_LIBRARIES})

# endif()

