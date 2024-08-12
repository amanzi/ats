# -*- mode: cmake -*-

#
#  ATS
#    density constitutive relations
#

add_library(lake_thermo_relations_density
  density_evaluator.cc
             )

install(TARGETS lake_thermo_relations_density DESTINATION lib)
