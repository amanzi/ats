# -*- mode: cmake -*-

#
#  ATS
#    Lake Thermo Heat Flux Evaluator class
#

include_directories(${ATS_SOURCE_DIR}/src/factory)

add_library(lake_thermo_relations_heat_flux_bc
            heat_flux_bc_evaluator.cc
            )

install(TARGETS lake_thermo_relations_heat_flux_bc DESTINATION lib)
