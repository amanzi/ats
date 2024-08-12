# -*- mode: cmake -*-

#
#  ATS
#    surface temperature constitutive relations
#

add_library(lake_thermo_relations_surface_temperature
  lake_surface_temperature_evaluator.cc
             )

install(TARGETS lake_thermo_relations_surface_temperature DESTINATION lib)
