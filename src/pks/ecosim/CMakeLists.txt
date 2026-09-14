# -*- mode: cmake -*-
#Everything here depends on EcoSIM
#so we put everything in an if statement

# -- not necessary it's only added if ENABLE_ECOSIM is set
if(ENABLE_ECOSIM)
    add_subdirectory(constitutive_relations)
    add_subdirectory(data)

    # These things should be set by FindEcoSIM.cmake
    set(ECOSIM_INSTALL_PREFIX ${ECOSIM_DIR})
    set(ECOSIM_INCLUDE_DIRS ${ECOSIM_DIR}/include)
    set(ECOSIM_LIB_DIR ${ECOSIM_DIR}/lib)

    message("In ATS-EcoSIM CMakeLists:")
    message("ECOSIM_DIR:" ${ECOSIM_DIR})
    message("ECOSIM_INSTALL_PREFIX:" ${ECOSIM_INSTALL_PREFIX})
    message("ECOSIM_INCLUDE_DIRS: ${ECOSIM_INCLUDE_DIRS}")
    message("ECOSIM_LIB_DIR:" ${ECOSIM_LIB_DIR})
    
    message("NetCDF_INCLUDE_DIRS: ${NetCDF_INCLUDE_DIRS}")
    message("NetCDF_C_LIBRARIES: ${NetCDF_C_LIBRARIES}")

    include_directories(${ATS_SOURCE_DIR}/src/pks)
    include_directories(${ATS_SOURCE_DIR}/src/pks/ecosim)
    include_directories(${ATS_SOURCE_DIR}/src/pks/ecosim/data)
    #include_directories(${ATS_SOURCE_DIR}/src/pks/ecosim/constitutive_relations/bulk_density)
    include_directories(${ATS_SOURCE_DIR}/src/pks/ecosim/constitutive_relations/hydraulic_conductivity)
    #include_directories(${ATS_SOURCE_DIR}/src/pks/ecosim/constitutive_relations/matric_pressure)

    include_directories(${ECOSIM_INCLUDE_DIRS})

    message("EcoSIM inc dirs: ")
    message(eco_inc_file="${ECOSIM_INCLUDE_DIRS}")

    # ATS EcoSIM pk
    # For adding F90 files add the F90 file to the ecosim source files and
    # inc files. The inc file also needs a header

    #set(ats_ecosim_src_files
    #    EcoSIM_ATS_interface.cc
    #    BGCEngine.cc
    #    ecosim_wrappers.F90
    #    data/bgc_fortran_memory_mod.F90
    #)

    #testing using link libs:
    #file(GLOB ECOLIBS ${ECOSIM_LIB_LOCATION}/*.a)

    set(ats_ecosim_src_files
      EcoSIM_ATS_interface.cc
      BGCEngine.cc
      data/bgc_fortran_memory_mod.F90
    )

    set(ats_ecosim_inc_files
      EcoSIM_ATS_interface.hh
      BGCEngine.hh
      ecosim_interface.h
  )

    # This should be set by FindEcoSIM.cmake
    # - not portable as *.a or *.dylib, hack for now *
    file(GLOB ECOSIM_LIBRARIES
      ${ECOSIM_LIB_DIR}/*
  )

    # not sure why this list is so long, and gfortran shouldn't be necessary
    # maybe we haven't identified the language correctly?
    set(ats_ecosim_link_libs
      ${Teuchos_LIBRARIES}
      ${Epetra_LIBRARIES}
      ${ECOSIM_LIBRARIES}
      # These are consistent with ats_pks_link_libs
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
      # Not sure about these ones?
      chemistry_pk
      ats_pks
      ats_eos
      ats_operators
      ats_ecosim_data
      ats_ecosim_relations
      ${NetCDF_C_LIBRARIES}
   )

    message(STATUS "ats_ecosim_link_libs: ${ats_ecosim_link_libs}")

    message("ats_ecosim_inc_files = ${ats_ecosim_inc_files}")
    message("ats_ecosim_link_libs = ${ats_ecosim_link_libs}")

    add_amanzi_library(ats_ecosim
                       SOURCE ${ats_ecosim_src_files}
                       HEADERS ${ats_ecosim_inc_files}
                       LINK_LIBS ${ats_ecosim_link_libs}
                       )

    #================================================
    # register evaluators/factories/pks

    register_evaluator_with_factory(
  HEADERFILE EcoSIM_ATS_interface_reg.hh
  LISTNAME ATS_ECOSIM_REG
  )

    generate_evaluators_registration_header(
  HEADERFILE ats_ecosim_registration.hh
  LISTNAME   ATS_ECOSIM_REG
  INSTALL    True
  )
endif()
