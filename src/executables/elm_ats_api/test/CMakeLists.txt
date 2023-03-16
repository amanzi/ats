# -*- mode: cmake -*-
#
#  ELM_ATS_API
#    example/test drivers
#
# build both tests for now
# need to make configurable

include_directories(${ATS_SOURCE_DIR}/src/executables)
include_directories(${ATS_SOURCE_DIR}/src/executables/elm_ats_api)

add_amanzi_executable(elm_ats_cc_test
  SOURCE elm_ats_test.cc
  LINK_LIBS elm_ats ${fates_link_libs} ${tpl_link_libs} ${ats_link_libs} ${amanzi_link_libs}
  OUTPUT_NAME elm_ats_cc_test
  OUTPUT_DIRECTORY ${ATS_BINARY_DIR})

add_subdirectory(fortran)
