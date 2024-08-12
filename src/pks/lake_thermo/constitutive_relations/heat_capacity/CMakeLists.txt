# -*- mode: cmake -*-

#
#  ATS
#    Lake Thermo PK class
#

include_directories(${ATS_SOURCE_DIR}/src/factory)

add_library(lake_thermo_relations_heat_capacity
            lake_heat_capacity_evaluator.cc
            )

install(TARGETS lake_thermo_relations_heat_capacity DESTINATION lib)
