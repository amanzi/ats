# -*- mode: cmake -*-

#
#  ATS
#    Flow PK class
#

# ATS include directories

add_subdirectory(constitutive_relations)
#add_subdirectory(morphology_pk)

include_directories(${GEOCHEM_SOURCE_DIR})
include_directories(${CHEMPK_SOURCE_DIR})
include_directories(${MPC_TREE_SOURCE_DIR})
include_directories(${ATS_SOURCE_DIR}/pks)
include_directories(${ATS_SOURCE_DIR}/constitutive_relations/eos)
include_directories(${ATS_SOURCE_DIR}/pks/flow)
include_directories(${ATS_SOURCE_DIR}/pks/energy)
include_directories(${ATS_SOURCE_DIR}/pks/transport)
include_directories(${ATS_SOURCE_DIR}/constitutive_relations/generic_evaluators)
include_directories(${ATS_SOURCE_DIR}/pks/flow/constitutive_relations/wrm)
include_directories(${ATS_SOURCE_DIR}/pks/flow/constitutive_relations/porosity)
include_directories(${ATS_SOURCE_DIR}/operators/upwinding)
include_directories(${ATS_SOURCE_DIR}/operators/advection)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/constitutive_relations)

set(ats_mpc_src_files
  weak_mpc.cc
  DomainSetMPC.cc
  operator_split_mpc.cc
  weak_mpc_semi_coupled.cc
  weak_mpc_semi_coupled_deform.cc
  mpc_coupled_cells.cc
  pk_mpcsubcycled_ats.cc
  mpc_surface_subsurface_helpers.cc
  weak_mpc_semi_coupled_helper.cc
  mpc_weak_subgrid.cc
  mpc_delegate_ewc.cc
  mpc_reactivetransport_pk.cc
  mpc_flowreactivetransport_pk.cc
  mpc_coupled_reactive_transport_pk.cc
  mpc_delegate_ewc_subsurface.cc
  mpc_delegate_water.cc
  mpc_coupled_water.cc
  mpc_coupled_water_split_flux.cc
  mpc_coupled_transport.cc
  mpc_permafrost_split_flux.cc
  mpc_permafrost_split_flux_columns.cc
  mpc_permafrost_split_flux_columns_subcycled.cc
  mpc_subsurface.cc
  mpc_surface.cc
  mpc_permafrost.cc
  mpc_morphology_pk.cc
  biomass_evaluator.cc
  )

set(ats_mpc_inc_files
  mpc.hh
  weak_mpc.hh
  strong_mpc.hh
  DomainSetMPC.hh
  operator_split_mpc.hh
  weak_mpc_semi_coupled.hh
  weak_mpc_semi_coupled_deform.hh
  mpc_coupled_cells.hh
  pk_mpcsubcycled_ats.hh
  mpc_surface_subsurface_helpers.hh
  weak_mpc_semi_coupled_helper.hh
  mpc_weak_subgrid.hh
  mpc_delegate_ewc.hh
  mpc_reactivetransport_pk.hh
  mpc_flowreactivetransport_pk.hh
  mpc_coupled_reactive_transport_pk.hh
  mpc_delegate_ewc_subsurface.hh
  mpc_delegate_water.hh
  mpc_coupled_water.hh
  mpc_coupled_transport.hh
  mpc_coupled_water_split_flux.hh
  mpc_permafrost_split_flux.hh
  mpc_permafrost_split_flux_columns.hh
  mpc_permafrost_split_flux_columns_subcycled.hh
  mpc_subsurface.hh
  mpc_surface.hh
  mpc_permafrost.hh
  biomass_evaluator.hh
  )

set(ats_mpc_link_libs
  ${Teuchos_LIBRARIES}
  ${Epetra_LIBRARIES}
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
  ats_operators
  ats_eos
  ats_generic_evals
  ats_pks
  ats_transport
  ats_flow
  ats_mpc_relations
  )

add_amanzi_library(ats_mpc
                   SOURCE ${ats_mpc_src_files}
                   HEADERS ${ats_mpc_inc_files}
		   LINK_LIBS ${ats_mpc_link_libs})

# register factories
register_evaluator_with_factory(
  HEADERFILE weak_mpc_reg.hh
  LISTNAME   ATS_MPC_REG
  )

register_evaluator_with_factory(
  HEADERFILE strong_mpc_reg.hh
  LISTNAME   ATS_MPC_REG
  )

register_evaluator_with_factory(
  HEADERFILE DomainSetMPC_reg.hh
  LISTNAME   ATS_MPC_REG
  )

register_evaluator_with_factory(
  HEADERFILE mpc_predictor_corrector_reg.hh
  LISTNAME   ATS_MPC_REG
  )

register_evaluator_with_factory(
  HEADERFILE operator_split_mpc_reg.hh
  LISTNAME   ATS_MPC_REG
  )

register_evaluator_with_factory(
  HEADERFILE pk_mpcsubcycled_ats_reg.hh
  LISTNAME   ATS_MPC_REG
  )

register_evaluator_with_factory(
   HEADERFILE mpc_flowreactivetransport_pk_reg.hh
   LISTNAME   ATS_MPC_REG
   )

register_evaluator_with_factory(
   HEADERFILE mpc_reactivetransport_pk_reg.hh
   LISTNAME   ATS_MPC_REG
   )
 
register_evaluator_with_factory(
  HEADERFILE mpc_coupled_reactive_transport_pk_reg.hh
  LISTNAME   ATS_MPC_REG
  )

register_evaluator_with_factory(
   HEADERFILE mpc_morphology_reg.hh
   LISTNAME   ATS_MPC_REG
   )
    
register_evaluator_with_factory(
  HEADERFILE mpc_subsurface_reg.hh
  LISTNAME   ATS_MPC_REG
  )

register_evaluator_with_factory(
  HEADERFILE mpc_surface_reg.hh
  LISTNAME   ATS_MPC_REG
  )

register_evaluator_with_factory(
  HEADERFILE mpc_coupled_cells_reg.hh
  LISTNAME   ATS_MPC_REG
  )

register_evaluator_with_factory(
  HEADERFILE mpc_permafrost_reg.hh
  LISTNAME   ATS_MPC_REG
  )

register_evaluator_with_factory(
  HEADERFILE mpc_coupled_water_reg.hh
  LISTNAME   ATS_MPC_REG
  )

register_evaluator_with_factory(
  HEADERFILE mpc_coupled_water_split_flux_reg.hh
  LISTNAME   ATS_MPC_REG
  )

register_evaluator_with_factory(
  HEADERFILE mpc_permafrost_split_flux_reg.hh
  LISTNAME   ATS_MPC_REG
  )

register_evaluator_with_factory(
  HEADERFILE mpc_permafrost_split_flux_columns_reg.hh
  LISTNAME   ATS_MPC_REG
  )

register_evaluator_with_factory(
  HEADERFILE mpc_permafrost_split_flux_columns_subcycled_reg.hh
  LISTNAME   ATS_MPC_REG
  )

register_evaluator_with_factory(
    HEADERFILE mpc_coupled_transport_reg.hh
    LISTNAME   ATS_MPC_REG
   )


register_evaluator_with_factory(
  HEADERFILE weak_mpc_semi_coupled_reg.hh
  LISTNAME   ATS_MPC_REG
  )

register_evaluator_with_factory(
  HEADERFILE mpc_weak_subgrid_reg.hh
  LISTNAME   ATS_MPC_REG
  )

register_evaluator_with_factory(
  HEADERFILE weak_mpc_semi_coupled_deform_reg.hh
  LISTNAME   ATS_MPC_REG
  )

register_evaluator_with_factory(
  HEADERFILE biomass_evaluator_reg.hh
  LISTNAME ATS_MPC_REG
  )

generate_evaluators_registration_header(
  HEADERFILE ats_mpc_registration.hh
  LISTNAME   ATS_MPC_REG
  INSTALL    True
)


