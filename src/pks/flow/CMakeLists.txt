# -*- mode: cmake -*-

#
#  ATS
#    Flow PK class
#

# ATS include directories
include_directories(${ATS_SOURCE_DIR}/src/pks)
include_directories(${ATS_SOURCE_DIR}/src/operators/advection)
include_directories(${ATS_SOURCE_DIR}/src/operators/upwinding)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/water_content)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/porosity)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/wrm)
# include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/overland_conductivity)
# include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/elevation)
# include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/sources)
# include_directories(${ATS_SOURCE_DIR}/src/constitutive_relations/eos)


set(ats_flow_src_files
  constitutive_relations/wrm/wrm_van_genuchten.cc
  richards_pk.cc
  richards_ti.cc
  richards_physics.cc
  richards_steadystate.cc
)

                 
  
  # predictor_delegate_bc_flux.cc
  # permafrost_pk.cc
  # interfrost.cc
  # overland_pressure_pk.cc
  # overland_pressure_physics.cc
  # overland_pressure_ti.cc
  #  overland_pk.cc
  #  overland_physics.cc
  #  overland_ti.cc
  # icy_overland.cc
  #  snow_distribution_pk.cc
  #  snow_distribution_physics.cc
  #  snow_distribution_ti.cc
  # preferential_pk.cc

set(ats_flow_inc_files
  constitutive_relations/water_content/richards_water_content_model.hh
  constitutive_relations/porosity/compressible_porosity_linear_model.hh
  constitutive_relations/wrm/wrm_van_genuchten.hh
  constitutive_relations/wrm/wrm_model.hh
  constitutive_relations/wrm/relative_permeability_model.hh

  richards.hh
  richards_steadystate.hh
)

  # flow_bc_factory.hh
  # predictor_delegate_bc_flux.hh
  # permafrost.hh
  # interfrost.hh
  # overland_pressure.hh
  # overland.hh
  # icy_overland.hh
  # snow_distribution.hh
  # preferential.hh


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
  # # ats_eos
  ats_pks
  )

add_amanzi_library(ats_flow
                   SOURCE ${ats_flow_src_files}
                   HEADERS ${ats_flow_inc_files}
		   LINK_LIBS ${ats_flow_link_libs})


register_abs_evaluator_with_factory(HEADERFILE "constitutive_relations/water_content/flow_relations_water_content_reg.hh"
  LISTNAME ATS_FLOW_RELATIONS_REG)
register_abs_evaluator_with_factory(HEADERFILE "constitutive_relations/porosity/flow_relations_porosity_reg.hh"
  LISTNAME ATS_FLOW_RELATIONS_REG)
register_abs_evaluator_with_factory(HEADERFILE "constitutive_relations/wrm/flow_relations_wrm_reg.hh"
  LISTNAME ATS_FLOW_RELATIONS_REG)

generate_evaluators_registration_header(
    HEADERFILE ats_flow_relations_registration.hh
    LISTNAME   ATS_FLOW_RELATIONS_REG
    INSTALL    True
    )

                 
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

# register_evaluator_with_factory(
#  HEADERFILE preferential_pk_reg.hh
#  LISTNAME   ATS_FLOW_PKS_REG
#  )

generate_evaluators_registration_header(
 HEADERFILE ats_flow_pks_registration.hh
 LISTNAME   ATS_FLOW_PKS_REG
 INSTALL    True
 )



 if (BUILD_TESTS)

    include_directories(${UnitTest_INCLUDE_DIRS})
    set(ats_flow_test_link_libs state mesh_factory mesh data_structures mesh_functions functions
                              ${UnitTest_LIBRARIES})

    add_amanzi_test(ats_flow_relations ats_flow_relations
     KIND int
     SOURCE test/Main.cc test/ats_flow_relations.cc
     LINK_LIBS ats_flow ${ats_flow_test_link_libs})

 endif()

