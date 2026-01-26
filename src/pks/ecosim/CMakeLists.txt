# -*- mode: cmake -*-
#Everything here depends on EcoSIM
#so we put everything in an if statement

if(ENABLE_ECOSIM)
    add_subdirectory(constitutive_relations)
    add_subdirectory(data)

    get_property(AMANZI_TPLS_DIR GLOBAL PROPERTY AMANZI_TPLS_DIR)

    set(ECOSIM_INSTALL_PREFIX ${ECOSIM_DIR}/ecosim)
    set(ECOSIM_LIB_LOCATION ${ECOSIM_DIR}/ecosim/local/lib)
    set(ECOSIM_BUILD_PREFIX ${ECOSIM_DIR}/ecosim/build)
    #set(ECOSIM_CMAKE_BINARY_DIR ${ECOSIM_DIR}/ecosim/build/Linux-x86_64-static-not-set-mpicc-Release)
    set(NETCDF_LIB ${ECOSIM_DIR}/lib)

    message("In ATS-EcoSIM CMakeLists:")
    message("ECOSIM_DIR:" ${ECOSIM_DIR})
    message("ECOSIM_INSTALL_PREFIX:" ${ECOSIM_INSTALL_PREFIX})
    message("ECOSIM_LIB_LOCATION:" ${ECOSIM_LIB_LOCATION})
    message("ECOSIM_BUILD_PREFIX:" ${ECOSIM_BUILD_PREFIX})

    include_directories(${ATS_SOURCE_DIR}/src/pks)
    include_directories(${ATS_SOURCE_DIR}/src/pks/ecosim)
    include_directories(${ATS_SOURCE_DIR}/src/pks/ecosim/data)
    #include_directories(${ATS_SOURCE_DIR}/src/pks/ecosim/constitutive_relations/bulk_density)
    include_directories(${ATS_SOURCE_DIR}/src/pks/ecosim/constitutive_relations/hydraulic_conductivity)
    #include_directories(${ATS_SOURCE_DIR}/src/pks/ecosim/constitutive_relations/matric_pressure)

    message("At include_directories")
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
    include_directories(${ECOSIM_BUILD_PREFIX}/f90src/HydroTherm/PhysData/)
    include_directories(${ECOSIM_BUILD_PREFIX}/f90src/HydroTherm/SnowPhys/)
    include_directories(${ECOSIM_BUILD_PREFIX}/f90src/Ecosim_mods/)
    include_directories(${ECOSIM_BUILD_PREFIX}/f90src/Prescribed_pheno/)
    include_directories(${ECOSIM_BUILD_PREFIX}/f90src/Plant_bgc/)
    include_directories(${ECOSIM_BUILD_PREFIX}/f90src/APIs/)
    include_directories(${ECOSIM_BUILD_PREFIX}/f90src/APIData/)

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
    ${ECOSIM_INSTALL_PREFIX}/f90src/ATSUtils/BGC_containers.F90
    ${ECOSIM_INSTALL_PREFIX}/f90src/ATSUtils/ATSCPLMod.F90
    ${ECOSIM_INSTALL_PREFIX}/f90src/ATSUtils/ATSEcoSIMInitMod.F90
    ${ECOSIM_INSTALL_PREFIX}/f90src/ATSUtils/ATSEcoSIMAdvanceMod.F90
    ${ECOSIM_INSTALL_PREFIX}/f90src/ATSUtils/ATSUtilsMod.F90
    ${ECOSIM_INSTALL_PREFIX}/f90src/ATSUtils/SharedDataMod.F90
    ${ECOSIM_INSTALL_PREFIX}/f90src/ATSUtils/c_f_interface_module.F90
    ${ECOSIM_INSTALL_PREFIX}/f90src/ATSUtils/ecosim_wrappers.F90
    )

    set(ats_ecosim_inc_files
    EcoSIM_ATS_interface.hh
    BGCEngine.hh
    ecosim_interface.h
  )

    file(GLOB ECOSIM_LIBRARIES
  ${ECOSIM_LIB_LOCATION}/*.a
  )


    find_package(NetCDF REQUIRED)

    set(ats_ecosim_link_libs
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
  ats_ecosim_data
  ats_ecosim_relations
  gfortran
  )

    message(STATUS "ats_ecosim_link_libs: ${ats_ecosim_link_libs}")

    message(inc_files="${ats_ecosim_inc_files}")

    add_amanzi_library(ats_ecosim
                    SOURCE ${ats_ecosim_src_files}
                    HEADERS ${ats_ecosim_inc_files}
        LINK_LIBS ${ats_ecosim_link_libs})

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
