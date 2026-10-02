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

# RPF - initial version attempted to use the same F90 file only
# in one place for ATS and E3SM. But, locating it only in the E3SM submodule
# seemed a bit too complicated for now and might be a better strategy if 
# one becomes a submodule of the other at any point
set(elm_water_balance_kernel_src "${CMAKE_CURRENT_SOURCE_DIR}/WaterBalanceKernel.F90")
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

if (BUILD_TESTS)

  # Unit test for the shared ELM water-balance kernel used by ATS.  It compiles
  # the *same* dependency-free WaterBalanceKernel.F90 the elm_ats library uses
  # (via elm_water_balance_kernel_src) and calls the bind(C) entry point
  # elm_ats_water_balance_error_c directly -- no ATS mesh/state/XML/Kokkos/MPI
  # needed.  This both checks ELM's balance formula and guards the C ABI symbol
  # that ELM_ATSDriver::checkELMWaterBalance_ relies on.
  #
  # The test is pure C++ with NO third-party dependency (no UnitTest++): it uses
  # a tiny self-contained harness and signals failure via a nonzero exit code,
  # which is all ctest needs.  Only the header dir is added so the test can find
  # elm_balance_interface_private.hh.
  include_directories(${CMAKE_CURRENT_SOURCE_DIR})

  add_amanzi_test(elm_ats_water_balance_kernel elm_ats_water_balance_kernel
    KIND unit
    SOURCE test/water_balance_kernel_test.cc
           ${elm_water_balance_kernel_src})

  # On Apple the mixed C++/Fortran link may need the same undefined-symbol
  # relaxation the shared library uses for the Fortran runtime.
  if (APPLE)
    set_target_properties(elm_ats_water_balance_kernel
      PROPERTIES LINK_FLAGS "-Wl,-undefined,dynamic_lookup")
  endif()

  # Unit test for the TimeAdvancer reject-and-retry path taken when the step
  # validity check (the ELM water-balance gate) fails.  Uses a stub PK and the
  # real kernel with a synthetic leak, so no mesh or input file is needed.
  include_directories(${ATS_SOURCE_DIR}/src/pks)

  add_amanzi_test(elm_ats_time_advancer_validity_retry elm_ats_time_advancer_validity_retry
    KIND unit
    SOURCE test/time_advancer_validity_retry_test.cc
           ${elm_water_balance_kernel_src}
    LINK_LIBS ${tpl_link_libs} ${ats_link_libs} ${amanzi_link_libs})

  if (APPLE)
    set_target_properties(elm_ats_time_advancer_validity_retry
      PROPERTIES LINK_FLAGS "-Wl,-undefined,dynamic_lookup")
  endif()

endif()
