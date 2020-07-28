# -*- mode: cmake -*-

#
#  ATS
#    WRM for flow
#

add_subdirectory(models)

include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/wrm/models)

add_library(flow_relations_wrm
  wrm_partition.cc
  wrm_permafrost_evaluator.cc
  wrm_evaluator.cc
  rel_perm_evaluator.cc
  pc_liquid_evaluator.cc
  pc_ice_evaluator.cc
  suction_head_evaluator.cc
  )

install(TARGETS flow_relations_wrm DESTINATION lib)
