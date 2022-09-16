# -*- mode: cmake -*-

set(fort_src_files
  ats_mod.f90
  test.f90
  )

add_amanzi_executable(fort_driver
  SOURCE ${fort_src_files}
  LINK_LIBS elm_ats_static
  OUTPUT_NAME elm_ats_f90_test
  OUTPUT_DIRECTORY ${ATS_BINARY_DIR})

