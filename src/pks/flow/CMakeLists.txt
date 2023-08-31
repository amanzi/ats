# -*- mode: cmake -*-

#
#  ATS
#    Flow PK class
#

add_subdirectory(constitutive_relations)

# ATS include directories
include_directories(${ATS_SOURCE_DIR}/src/pks)
include_directories(${ATS_SOURCE_DIR}/src/operators/advection)
include_directories(${ATS_SOURCE_DIR}/src/operators/upwinding)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/water_content)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/wrm)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/overland_conductivity)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/elevation)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/sources)
include_directories(${ATS_SOURCE_DIR}/src/constitutive_relations/eos)


set(ats_flow_src_files
  predictor_delegate_bc_flux.cc
  richards_pk.cc
  richards_ti.cc
  richards_physics.cc
  richards_steadystate.cc
  permafrost_pk.cc
  interfrost.cc
  overland_pressure_pk.cc
  overland_pressure_physics.cc
  overland_pressure_ti.cc
  #  overland_pk.cc
  #  overland_physics.cc
  #  overland_ti.cc
  icy_overland.cc
  snow_distribution_pk.cc
  snow_distribution_physics.cc
  snow_distribution_ti.cc
  preferential_pk.cc
  )

set(ats_flow_inc_files
  flow_bc_factory.hh
  predictor_delegate_bc_flux.hh
  richards.hh
  richards_steadystate.hh
  permafrost.hh
  interfrost.hh
  overland_pressure.hh
  overland.hh
  icy_overland.hh
  snow_distribution.hh
  preferential.hh
  )

set(ats_flow_link_libs
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
  ats_operators
  ats_eos
  ats_pks
  ats_flow_relations
  )

add_amanzi_library(ats_flow
                   SOURCE ${ats_flow_src_files}
                   HEADERS ${ats_flow_inc_files}
		   LINK_LIBS ${ats_flow_link_libs})


#
# generate registration files
#
register_evaluator_with_factory(
  HEADERFILE richards_pk_reg.hh
  LISTNAME   ATS_FLOW_PKS_REG
  )

register_evaluator_with_factory(
  HEADERFILE richards_steadystate_reg.hh
  LISTNAME   ATS_FLOW_PKS_REG
  )

register_evaluator_with_factory(
  HEADERFILE permafrost_pk_reg.hh
  LISTNAME   ATS_FLOW_PKS_REG
  )

register_evaluator_with_factory(
  HEADERFILE interfrost_pk_reg.hh
  LISTNAME   ATS_FLOW_PKS_REG
  )

# register_evaluator_with_factory(
#   HEADERFILE overland_pk_reg.hh
#   LISTNAME   ATS_FLOW_PKS_REG
#   )

register_evaluator_with_factory(
  HEADERFILE overland_pressure_pk_reg.hh
  LISTNAME   ATS_FLOW_PKS_REG
  )

register_evaluator_with_factory(
  HEADERFILE icy_overland_reg.hh
  LISTNAME   ATS_FLOW_PKS_REG
  )

 register_evaluator_with_factory(
   HEADERFILE snow_distribution_reg.hh
   LISTNAME   ATS_FLOW_PKS_REG
   )

register_evaluator_with_factory(
 HEADERFILE preferential_pk_reg.hh
 LISTNAME   ATS_FLOW_PKS_REG
 )

generate_evaluators_registration_header(
 HEADERFILE ats_flow_pks_registration.hh
 LISTNAME   ATS_FLOW_PKS_REG
 INSTALL    True
 )





