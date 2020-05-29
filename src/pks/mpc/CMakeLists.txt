# -*- mode: cmake -*-

#
#  ATS
#    Flow PK class
#

# ATS include directories
include_directories(${ATS_SOURCE_DIR}/pks/flow)

set(ats_mpc_src_files
  )

set(ats_mpc_inc_files
  mpc_coupled_water.hh
  )

# set(ats_mpc_inc_files
#   mpc.hh
#   weak_mpc.hh
#   strong_mpc.hh
#   operator_split_mpc.hh
#   weak_mpc_semi_coupled.hh
#   weak_mpc_semi_coupled_deform.hh
#   mpc_coupled_cells.hh
#   pk_mpcsubcycled_ats.hh
#   mpc_surface_subsurface_helpers.hh
#   weak_mpc_semi_coupled_helper.hh
#   mpc_weak_subgrid.hh
#   mpc_delegate_ewc.hh
#   mpc_reactivetransport_pk.hh
#   mpc_flowreactivetransport_pk.hh
#   mpc_delegate_ewc_subsurface.hh
#   mpc_delegate_water.hh
#   mpc_coupled_water.hh
#   mpc_coupled_transport.hh
#   mpc_coupled_water_split_flux.hh
#   mpc_permafrost_split_flux.hh
#   mpc_permafrost_split_flux_columns.hh
#   mpc_permafrost_split_flux_columns_subcycled.hh
#   mpc_subsurface.hh
#   mpc_surface.hh
#   mpc_permafrost.hh
#   biomass_evaluator.hh
#   )

set(ats_mpc_link_libs
  ${Teuchos_LIBRARIES}
  ${Tpetra_LIBRARIES}
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
  )

# add_amanzi_library(ats_mpc
#                    SOURCE ${ats_mpc_src_files}
#                    HEADERS ${ats_mpc_inc_files}
# 		   LINK_LIBS ${ats_mpc_link_libs})

# register factories
register_evaluator_with_factory(
  HEADERFILE mpc_coupled_water_reg.hh
  LISTNAME   ATS_MPC_REG
  )

# register_evaluator_with_factory(
#   HEADERFILE operator_split_mpc_reg.hh
#   LISTNAME   ATS_MPC_REG
#   )

# register_evaluator_with_factory(
#   HEADERFILE pk_mpcsubcycled_ats_reg.hh
#   LISTNAME   ATS_MPC_REG
#   )

# register_evaluator_with_factory(
#    HEADERFILE mpc_flowreactivetransport_pk_reg.hh
#    LISTNAME   ATS_MPC_REG
#    )

# register_evaluator_with_factory(
#    HEADERFILE mpc_reactivetransport_pk_reg.hh
#    LISTNAME   ATS_MPC_REG
#    )

# register_evaluator_with_factory(
#    HEADERFILE mpc_morphology_reg.hh
#    LISTNAME   ATS_MPC_REG
#    )
    
# register_evaluator_with_factory(
#   HEADERFILE mpc_subsurface_reg.hh
#   LISTNAME   ATS_MPC_REG
#   )

# register_evaluator_with_factory(
#   HEADERFILE mpc_surface_reg.hh
#   LISTNAME   ATS_MPC_REG
#   )

# register_evaluator_with_factory(
#   HEADERFILE mpc_coupled_cells_reg.hh
#   LISTNAME   ATS_MPC_REG
#   )

# register_evaluator_with_factory(
#   HEADERFILE mpc_permafrost_reg.hh
#   LISTNAME   ATS_MPC_REG
#   )

# register_evaluator_with_factory(
#   HEADERFILE mpc_coupled_water_reg.hh
#   LISTNAME   ATS_MPC_REG
#   )

# register_evaluator_with_factory(
#   HEADERFILE mpc_coupled_water_split_flux_reg.hh
#   LISTNAME   ATS_MPC_REG
#   )

# register_evaluator_with_factory(
#   HEADERFILE mpc_permafrost_split_flux_reg.hh
#   LISTNAME   ATS_MPC_REG
#   )

# register_evaluator_with_factory(
#   HEADERFILE mpc_permafrost_split_flux_columns_reg.hh
#   LISTNAME   ATS_MPC_REG
#   )

# register_evaluator_with_factory(
#   HEADERFILE mpc_permafrost_split_flux_columns_subcycled_reg.hh
#   LISTNAME   ATS_MPC_REG
#   )

# register_evaluator_with_factory(
#     HEADERFILE mpc_coupled_transport_reg.hh
#     LISTNAME   ATS_MPC_REG
#    )


# register_evaluator_with_factory(
#   HEADERFILE weak_mpc_semi_coupled_reg.hh
#   LISTNAME   ATS_MPC_REG
#   )

# register_evaluator_with_factory(
#   HEADERFILE mpc_weak_subgrid_reg.hh
#   LISTNAME   ATS_MPC_REG
#   )

# register_evaluator_with_factory(
#   HEADERFILE weak_mpc_semi_coupled_deform_reg.hh
#   LISTNAME   ATS_MPC_REG
#   )

# register_evaluator_with_factory(
#   HEADERFILE biomass_evaluator_reg.hh
#   LISTNAME ATS_MPC_REG
#   )

# generate_evaluators_registration_header(
#   HEADERFILE ats_mpc_registration.hh
#   LISTNAME   ATS_MPC_REG
#   INSTALL    True
# )




