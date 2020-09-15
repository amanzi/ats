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

#project(EXECUTABLE)

#
#  ATS
#    Executable
#
include_directories(${GEOCHEM_SOURCE_DIR})
include_directories(${MESH_FACTORY_SOURCE_DIR})
include_directories(${MESH_LOGICAL_SOURCE_DIR})
include_directories(${MSTK_SOURCE_DIR})
include_directories(${CHEMPK_SOURCE_DIR})
include_directories(${MPC_TREE_SOURCE_DIR})
include_directories(${DBG_SOURCE_DIR})
include_directories(${TRANSPORT_SOURCE_DIR})
get_property(CHEM_INCLUDES_DIR GLOBAL PROPERTY CHEM_INCLUDES_DIR)
include_directories(${CHEM_INCLUDES_DIR})

#include_directories(${ATS_SOURCE_DIR}/constitutive_relations)
#include_directories(${ATS_SOURCE_DIR}/constitutive_relations/eos)
#include_directories(${ATS_SOURCE_DIR}/constitutive_relations/surface_subsurface_fluxes)
#include_directories(${ATS_SOURCE_DIR}/constitutive_relations/generic_evaluators)
include_directories(${ATS_SOURCE_DIR}/pks)
#include_directories(${ATS_SOURCE_DIR}/pks/mpc)
#include_directories(${ATS_SOURCE_DIR}/pks/energy)
#include_directories(${ATS_SOURCE_DIR}/pks/flow)
#include_directories(${ATS_SOURCE_DIR}/pks/deformation)
#include_directories(${ATS_SOURCE_DIR}/pks/transport)
include_directories(${ATS_SOURCE_DIR}/operators/upwinding)
include_directories(${ATS_SOURCE_DIR}/operators/advection)
include_directories(${ATS_SOURCE_DIR}/operators/deformation)

if (ENABLE_FATES)
  link_directories(/nh/u/dasvyat/Coastal/ats-fates-new/install/lib/)
  link_directories(/nh/u/dasvyat/Coastal/ParallelIO/install/lib/)
  link_directories(/nh/u/dasvyat/Coastal/pnetcdf-1.11.2/install/lib)
  set (VEG_LIB pk_veg)
endif()


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
  coordinator.cc
  ats_mesh_factory.cc
  simulation_driver.cc
  main.cc
  )

set(ats_inc_files
  coordinator.hh
  ats_mesh_factory.hh
  simulation_driver.hh
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
  mesh_audit
  mesh_functions
  functions
  geometry
  mesh_factory
  output
  mstk_mesh
  mesh_logical
  chemistry_pk
  mpc_tree
  transport
  )
  
set(ats_link_libs
  ats_operators
  ats_generic_evals
  ats_surf_subsurf
  ats_eos
  ats_pks
  ats_transport
  ats_sed_transport  
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

if (ENABLE_FATES)
  set(fates_link_libs -lfates -lparteh -lshared  -lpiof -lpioc -lpnetcdf -lgptl)
#  link_directories(ats /nh/u/dasvyat/Coastal/ats-fates-new/install/lib/)
endif()


add_amanzi_executable(ats
  SOURCE ${ats_src_files}
  LINK_LIBS ${fates_link_libs} ${tpl_link_libs} ${ats_link_libs} ${amanzi_link_libs} 
  OUTPUT_NAME ats
  OUTPUT_DIRECTORY ${ATS_BINARY_DIR})
