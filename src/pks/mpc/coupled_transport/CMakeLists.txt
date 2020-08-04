# -*- mode: cmake -*-

#
#  ATS
#    Flow PK class
#

# ATS include directories


# include_directories(${ATS_SOURCE_DIR}/src/data_structures)
# include_directories(${ATS_SOURCE_DIR}/src/state)
# include_directories(${ATS_SOURCE_DIR}/src/constitutive_relations/eos)
# include_directories(${ATS_SOURCE_DIR}/src/constitutive_relations/ewc)
# include_directories(${ATS_SOURCE_DIR}/src/constitutive_relations/generic_evaluators)
# include_directories(${ATS_SOURCE_DIR}/src/pks)
# include_directories(${ATS_SOURCE_DIR}/src/pks/mpc)
# include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations)
# include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/porosity)
# include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/elevation)
# include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/wrm)
# include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/overland_conductivity)
# include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/water_content)
# include_directories(${ATS_SOURCE_DIR}/src/pks/flow/richards)
# include_directories(${ATS_SOURCE_DIR}/src/pks/flow/permafrost)
# include_directories(${ATS_SOURCE_DIR}/src/pks/flow/overland)
# include_directories(${ATS_SOURCE_DIR}/src/pks/flow/overland_pressure)
include_directories(${ATS_SOURCE_DIR}/src/pks/transport/transport_amanzi)


add_library(pk_coupled_transport mpc_coupled_transport.cc)

install(TARGETS pk_coupled_transport  DESTINATION lib)


include_evaluators_directories(LISTNAME ATS_POROSITY_EVAL_INCLUDES)
include_evaluators_directories(LISTNAME CONSTITUTIVE_RELATIONS_EOS_EVALUATORS_INCLUDES)
include_evaluators_directories(LISTNAME CONSTITUTIVE_RELATIONS_SURFACE_SUBSURFACE_FLUXES_EVALUATORS_INCLUDES)
include_evaluators_directories(LISTNAME CONSTITUTIVE_RELATIONS_GENERIC_EVALUATORS_INCLUDES)
include_evaluators_directories(LISTNAME FLOW_CONSTITUTIVE_RELATIONS_POROSITY_EVALUATORS_INCLUDES)
include_evaluators_directories(LISTNAME WRM_MODELS_EVALUATORS_FACTORIES_INCLUDES)
include_evaluators_directories(LISTNAME WMR_EVALUATORS_FACTORIES_INCLUDES)
include_evaluators_directories(LISTNAME RICHARDS_EVAL_FAC_INCLUDES)
include_evaluators_directories(LISTNAME MPC_PK_FAC_INCLUDES)


