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
#include_directories(${MPC_TREE_SOURCE_DIR})
include_directories(${TRANSPORT_SOURCE_DIR})
include_directories(${ATS_SOURCE_DIR}/src/pks)
include_directories(${ATS_SOURCE_DIR}/src/constitutive_relations/eos)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow)
include_directories(${ATS_SOURCE_DIR}/src/pks/energy)
include_directories(${ATS_SOURCE_DIR}/src/pks/transport)
include_directories(${ATS_SOURCE_DIR}/src/pks/surface_balance)
include_directories(${ATS_SOURCE_DIR}/src/constitutive_relations/generic_evaluators)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/wrm)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/porosity)
include_directories(${ATS_SOURCE_DIR}/src/pks/surface_balance/constitutive_relations/land_cover)
include_directories(${ATS_SOURCE_DIR}/src/operators/upwinding)
include_directories(${ATS_SOURCE_DIR}/src/operators/advection)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/constitutive_relations)

set(ats_mpc_src_files
  weak_mpc.cc
  mpc_subcycled.cc
  mpc_surface_subsurface_helpers.cc
  mpc_coupled_cells.cc
  mpc_delegate_ewc.cc
  mpc_delegate_ewc_subsurface.cc
  mpc_subsurface.cc
  mpc_delegate_ewc_surface.cc
  mpc_surface.cc
  mpc_delegate_water.cc
  mpc_coupled_water.cc
#  mpc_coupled_dualmedia_water.cc
  mpc_coupled_water_split_flux.cc
  mpc_permafrost.cc
  mpc_coupled_transport.cc
  mpc_reactivetransport.cc
  mpc_coupled_reactivetransport.cc
  mpc_flow_transport.cc
  mpc_weak_subdomain.cc
  mpc_coupled_water_split_flux.cc
  mpc_permafrost_split_flux.cc
  mpc_morphology_pk.cc
#  biomass_evaluator.cc
  )

set(ats_mpc_inc_files
  mpc.hh
  weak_mpc.hh
  strong_mpc.hh
  mpc_subcycled.hh
  mpc_surface_subsurface_helpers.hh
  mpc_coupled_cells.hh
  mpc_delegate_ewc.hh
  mpc_delegate_ewc_subsurface.hh
  mpc_subsurface.hh
  mpc_delegate_ewc_surface.hh
  mpc_surface.hh
  mpc_delegate_water.hh
  mpc_coupled_water.hh
  mpc_coupled_dualmedia_water.hh
  mpc_permafrost.hh
  mpc_coupled_transport.hh
  mpc_reactivetransport.hh
  mpc_coupled_reactivetransport.hh
  mpc_flow_transport.hh
  mpc_weak_subdomain.hh
  mpc_coupled_water_split_flux.hh
  mpc_permafrost_split_flux.hh
#  biomass_evaluator.hh
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
  ats_energy
  ats_flow
  ats_surface_balance
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
  HEADERFILE mpc_subcycled_reg.hh
  LISTNAME   ATS_MPC_REG
  )

register_evaluator_with_factory(
  HEADERFILE mpc_reactivetransport_reg.hh
  LISTNAME   ATS_MPC_REG
  )
 
register_evaluator_with_factory(
  HEADERFILE mpc_coupled_reactivetransport_reg.hh
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

# register_evaluator_with_factory(
#   HEADERFILE mpc_coupled_dualmedia_water_reg.hh
#   LISTNAME   ATS_MPC_REG
#   )


register_evaluator_with_factory(
  HEADERFILE mpc_coupled_water_split_flux_reg.hh
  LISTNAME   ATS_MPC_REG
  )

register_evaluator_with_factory(
  HEADERFILE mpc_permafrost_split_flux_reg.hh
  LISTNAME   ATS_MPC_REG
  )

register_evaluator_with_factory(
  HEADERFILE mpc_coupled_transport_reg.hh
  LISTNAME   ATS_MPC_REG
  )

register_evaluator_with_factory(
  HEADERFILE mpc_weak_subdomain_reg.hh
  LISTNAME   ATS_MPC_REG
  )

#register_evaluator_with_factory(
#  HEADERFILE biomass_evaluator_reg.hh
#  LISTNAME ATS_MPC_REG
#  )

generate_evaluators_registration_header(
  HEADERFILE ats_mpc_registration.hh
  LISTNAME   ATS_MPC_REG
  INSTALL    True
)
