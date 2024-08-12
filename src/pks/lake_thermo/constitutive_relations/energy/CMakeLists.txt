# -*- mode: cmake -*-

#
#  ATS
#    energy constitutive relations
#

add_library(lake_thermo_relations_energy
  lake_energy_evaluator.cc
             )

install(TARGETS lake_thermo_relations_energy DESTINATION lib)
