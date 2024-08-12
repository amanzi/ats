# -*- mode: cmake -*-

#
#  ATS
#    Lake Thermo PK class
#

include_directories(${ATS_SOURCE_DIR}/src/factory)

add_library(lake_thermo_relations_evaporation_rate
            lake_evaporation_rate_evaluator.cc
            )

install(TARGETS lake_thermo_relations_evaporation_rate DESTINATION lib)
