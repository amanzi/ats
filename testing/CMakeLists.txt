#
# This adds the ATS tests, which are in a separate repo, into the
# cmake build so that they build via ctests/make test
#
if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/ats-regression-tests/README.md")
  message("-- Adding ats-regression-tests")

  # Copy test subdirectory for out of source builds
  if (NOT ("${ATS_SOURCE_DIR}/src" STREQUAL "${ATS_BINARY_DIR}"))
    set(ATS_REGRESSION_TESTS_DIR "${ATS_BINARY_DIR}/testing/ats-regression-tests")
    file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/ats-regression-tests"
      DESTINATION "${ATS_BINARY_DIR}/testing")
  else()
    set(ATS_REGRESSION_TESTS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/ats-regression-tests")    
  endif()

  # figure out what tests are available
  if (ENABLE_Regression_Tests)
    set(ats_available_tests_string "")
    execute_process(COMMAND ${PYTHON_EXECUTABLE} ${ATS_REGRESSION_TESTS_DIR}/regression_tests.py --suites=testing --ats=${ATS_SOURCE_DIR} --list-tests .
      OUTPUT_VARIABLE ats_available_tests_string
      WORKING_DIRECTORY ${ATS_REGRESSION_TESTS_DIR})

    message("-- Found ATS tests:")
    string(STRIP ${ats_available_tests_string} ats_available_tests_string)
    string(REPLACE "\n" ";" ats_available_tests ${ats_available_tests_string})

    foreach(ats_test_and_dir_string IN LISTS ats_available_tests)
      # add a test that calls the python regression test suite
      message("   adding test: ${ats_test_and_dir_string}")
      string(REPLACE " " ";" ats_test_and_dir ${ats_test_and_dir_string})

      list(GET ats_test_and_dir 0 ats_test_dir)
      list(GET ats_test_and_dir 1 ats_test_name)
      set(ats_test_full_name "ats_regression_test-${ats_test_dir}-${ats_test_name}")
      
      add_test(NAME ${ats_test_full_name}
        COMMAND ${PYTHON_EXECUTABLE} ${ATS_REGRESSION_TESTS_DIR}/regression_tests.py -e ${ATS_BINARY_DIR}/src/executables/ats ${ats_test_dir} -t ${ats_test_name}
        WORKING_DIRECTORY ${ATS_REGRESSION_TESTS_DIR})

    endforeach()

  endif()
endif()

