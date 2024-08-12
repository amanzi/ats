# -*- mode: cmake -*-

#
#  ATS
#    enthalpy constitutive relations
#

add_library(lake_thermo_relations_enthalpy
  enthalpy_evaluator.cc
             )

install(TARGETS lake_thermo_relations_enthalpy DESTINATION lib)
