# -*- mode: cmake -*-

# # Need this define. Errors from MSTK include files 
# # about MPI_COMM_WORLD. --lpritch
# add_definitions("-DMSTK_HAVE_MPI")
# add_definitions("-DDISABLE_PHYSICS")

# if (WITH_MSTK_2_20rc1_OR_NEWER)
#   add_definitions("-DMSTK_2_20rc1_OR_NEWER")
# endif ()

# if (WITH_MSTK_2_21rc1_OR_NEWER)
#   add_definitions("-DMSTK_2_21rc1_OR_NEWER")
# endif ()

project(EXECUTABLE)

#
#  ATS
#    Executable
#
include_directories(${GEOCHEM_SOURCE_DIR})
include_directories(${MESH_FACTORY_SOURCE_DIR})
include_directories(${MESH_LOGICAL_SOURCE_DIR})
include_directories(${MESH_EXTRACTED_SOURCE_DIR})
include_directories(${MESH_MSTK_SOURCE_DIR})
include_directories(${CHEMPK_SOURCE_DIR})
include_directories(${MPC_TREE_SOURCE_DIR})
include_directories(${DBG_SOURCE_DIR})
include_directories(${TRANSPORT_SOURCE_DIR})
get_property(CHEM_INCLUDES_DIR GLOBAL PROPERTY CHEM_INCLUDES_DIR)
include_directories(${CHEM_INCLUDES_DIR})

include_directories(${ATS_SOURCE_DIR}/src/constitutive_relations)
include_directories(${ATS_SOURCE_DIR}/src/constitutive_relations/eos)
include_directories(${ATS_SOURCE_DIR}/src/constitutive_relations/surface_subsurface_fluxes)
include_directories(${ATS_SOURCE_DIR}/src/constitutive_relations/generic_evaluators)
include_directories(${ATS_SOURCE_DIR}/src/constitutive_relations/column_integrators)
include_directories(${ATS_SOURCE_DIR}/src/pks)
include_directories(${ATS_SOURCE_DIR}/src/pks/mpc)
include_directories(${ATS_SOURCE_DIR}/src/pks/energy)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow)
include_directories(${ATS_SOURCE_DIR}/src/pks/deformation)
include_directories(${ATS_SOURCE_DIR}/src/pks/transport)
include_directories(${ATS_SOURCE_DIR}/src/operators/upwinding)
include_directories(${ATS_SOURCE_DIR}/src/operators/advection)
include_directories(${ATS_SOURCE_DIR}/src/operators/deformation)

include_directories(${AMANZI_BINARY_DIR}) # required to pick up amanzi_version.hh
include_directories(${ATS_BINARY_DIR})

include_evaluators_directories(LISTNAME REGISTER_AMANZI_STATE_EVALUATORS_INCLUDES)
include_evaluators_directories(LISTNAME ATS_RELATIONS_REG_INCLUDES)
include_evaluators_directories(LISTNAME ATS_TRANSPORT_REG_INCLUDES)
include_evaluators_directories(LISTNAME ATS_ENERGY_PKS_REG_INCLUDES)
include_evaluators_directories(LISTNAME ATS_ENERGY_RELATIONS_REG_INCLUDES)
include_evaluators_directories(LISTNAME ATS_FLOW_PKS_REG_INCLUDES)
include_evaluators_directories(LISTNAME ATS_FLOW_RELATIONS_REG_INCLUDES)
include_evaluators_directories(LISTNAME ATS_DEFORMATION_REG_INCLUDES)
include_evaluators_directories(LISTNAME ATS_SURFACE_BALANCE_REG_INCLUDES)
include_evaluators_directories(LISTNAME ATS_BGC_REG_INCLUDES)
include_evaluators_directories(LISTNAME ATS_MPC_REG_INCLUDES)
include_evaluators_directories(LISTNAME SED_TRANSPORT_REG_INCLUDES)

set(ats_src_files
  ats_mesh_factory.cc
  coordinator.cc
  ats_driver.cc
  )

set(ats_inc_files
  ats_mesh_factory.hh
  coordinator.hh
  ats_driver.hh
  )

set(amanzi_link_libs
  operators
  pks
  geochemrxns
  geochembase
  geochemutil
  geochemsolvers
  state
  whetstone
  time_integration
  solvers
  dbg       
  data_structures
  mesh
  mesh_functions
  functions
  geometry
  mesh_factory
  output
  mesh_mstk
  mesh_logical
  chemistry_pk
  #mpc_tree
  transport
  )
  
set(ats_link_libs
  ats_operators
  ats_generic_evals
  ats_column_integrator
  ats_surf_subsurf
  ats_eos
  ats_pks
  ats_transport
  #ats_sed_transport  
  ats_energy
  ats_energy_relations
  ats_flow
  ats_flow_relations
  ats_deform
  ats_bgc
  ats_surface_balance
  ats_mpc
  ats_mpc_relations
  )


# note, we can be inclusive here, because if they aren't enabled,
# these won't be defined and will result in empty strings.
set(tpl_link_libs
  ${ALQUIMIA_LIBRARIES}
  ${PFLOTRAN_LIBRARIES}
  ${Teuchos_LIBRARIES}
  ${Epetra_LIBRARIES}
  ${Boost_LIBRARIES}
  ${PETSC_LIBRARIES}
  ${MSTK_LIBRARIES}
  ${SILO_LIBRARIES}
  ${HYPRE_LIBRARIES}
  ${HDF5_LIBRARIES}
  ${CLM_LIBRARIES}
  )

add_amanzi_library(ats_executable
  SOURCE ${ats_src_files}
  HEADERS ${ats_inc_files}
  LINK_LIBS ${amanzi_link_libs} ${tpl_link_libs})
if (APPLE AND BUILD_SHARED_LIBS)
  set_target_properties(ats_executable PROPERTIES LINK_FLAGS "-Wl,-undefined,dynamic_lookup")
endif()			  

if (BUILD_TESTS) 

  # Add UnitTest includes
  include_directories(${UnitTest_INCLUDE_DIRS})
  include_directories(${EXECUTABLE_SOURCE_DIR})

  # Copy test subdirectory for out of source builds
  if (NOT ("${EXECUTABLE_SOURCE_DIR}" STREQUAL "${EXECUTABLE_BINARY_DIR}"))
    file(GLOB DataFiles "${EXECUTABLE_SOURCE_DIR}/test/*.xml"
      "${EXECUTABLE_SOURCE_DIR}/test/*.exo"
      "${EXECUTABLE_SOURCE_DIR}/test/*.h5"
      "${EXECUTABLE_SOURCE_DIR}/test/*.bin")
    file(COPY ${DataFiles} DESTINATION ${EXECUTABLE_BINARY_DIR}/test/)
  endif()

  # test for ats mesh factory
  add_amanzi_test(executable_mesh_factory executable_mesh_factory
           KIND int
           SOURCE test/Main.cc test/executable_mesh_factory.cc
           LINK_LIBS ats_executable ${ats_link_libs} ${UnitTest_LIBRARIES} ${NOX_LIBRARIES} ${HDF5_LIBRARIES})
  add_amanzi_test(executable_mesh_factory_np2 executable_mesh_factory NPROCS 2 KIND uint)
  add_amanzi_test(executable_mesh_factory_np4 executable_mesh_factory NPROCS 2 KIND uint)

  # test for coupled water preconditioners
  add_amanzi_test(executable_coupled_water executable_coupled_water
    KIND int
    SOURCE test/Main.cc test/executable_coupled_water.cc
    LINK_LIBS ats_executable ${ats_link_libs} ${UnitTest_LIBRARIES} ${NOX_LIBRARIES} ${HDF5_LIBRARIES})

endif()

add_amanzi_executable(ats
  SOURCE main.cc
  LINK_LIBS ats_executable ${fates_link_libs} ${tpl_link_libs} ${ats_link_libs} ${amanzi_link_libs} 
  OUTPUT_NAME ats
  OUTPUT_DIRECTORY ${ATS_BINARY_DIR})
