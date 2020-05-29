# -*- mode: cmake -*-

#
#  ATS
#    Flow PK class
#

# ATS include directories
#include_directories(${ATS_SOURCE_DIR}/pks)
# include_directories(${ATS_SOURCE_DIR}/operators/advection)
# include_directories(${ATS_SOURCE_DIR}/operators/upwinding)
# include_directories(${ATS_SOURCE_DIR}/pks/flow/constitutive_relations/water_content)
# include_directories(${ATS_SOURCE_DIR}/pks/flow/constitutive_relations/wrm)
# include_directories(${ATS_SOURCE_DIR}/pks/flow/constitutive_relations/overland_conductivity)
# include_directories(${ATS_SOURCE_DIR}/pks/flow/constitutive_relations/elevation)
# include_directories(${ATS_SOURCE_DIR}/pks/flow/constitutive_relations/sources)
# include_directories(${ATS_SOURCE_DIR}/constitutive_relations/eos)


set(ats_flow_src_files
  )

set(ats_flow_inc_files
  models/CapillaryPressureLiqAtm.hh
  models/CompressiblePorosityLinear.hh
  models/ManningConductivity.hh
  models/SurfaceWaterContent.hh
  models/SurfaceWaterDepth.hh
  models/WRMVanGenuchten.hh
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
  )


# add_amanzi_library(ats_flow
#                    SOURCE ${ats_flow_src_files}
#                    HEADERS ${ats_flow_inc_files}
# 		   LINK_LIBS ${ats_flow_link_libs})


#
# generate registration files
#
register_evaluator_with_factory(
  HEADERFILE richards_reg.hh
  LISTNAME   ATS_REG
  )

register_evaluator_with_factory(
  HEADERFILE models/FlowModels_reg.hh
  LISTNAME  ATS_REG
  )


# register_evaluator_with_factory(
#   HEADERFILE richards_steadystate_reg.hh
#   LISTNAME   ATS_FLOW_PKS_REG
#   )

# register_evaluator_with_factory(
#   HEADERFILE permafrost_pk_reg.hh
#   LISTNAME   ATS_FLOW_PKS_REG
#   )

# register_evaluator_with_factory(
#   HEADERFILE interfrost_pk_reg.hh
#   LISTNAME   ATS_FLOW_PKS_REG
#   )

# register_evaluator_with_factory(
#   HEADERFILE overland_pk_reg.hh
#   LISTNAME   ATS_FLOW_PKS_REG
#   )

# register_evaluator_with_factory(
#   HEADERFILE overland_pressure_pk_reg.hh
#   LISTNAME   ATS_FLOW_PKS_REG
#   )

# register_evaluator_with_factory(
#   HEADERFILE icy_overland_reg.hh
#   LISTNAME   ATS_FLOW_PKS_REG
#   )

# register_evaluator_with_factory(
#   HEADERFILE snow_distribution_reg.hh
#   LISTNAME   ATS_FLOW_PKS_REG
#   )






