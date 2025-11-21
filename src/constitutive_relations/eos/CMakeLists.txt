# -*- mode: cmake -*-

#
#  ATS
#    Equations of state
#

set(ats_eos_src_files
  eos_factory.cc
  eos_evaluator.cc
  eos_constant.cc
  eos_linear.cc
  eos_ideal_gas.cc
  eos_water.cc
  eos_ice.cc
  eos_vapor_in_gas.cc
  eos_sw.cc
  viscosity_evaluator.cc
  viscosity_relation_factory.cc
  viscosity_constant.cc
  viscosity_water.cc
  molar_fraction_gas_evaluator.cc
  vapor_pressure_relation_factory.cc
  vapor_pressure_water.cc
  carbon_decomposition_rate_evaluator.cc
  transport_decay_rate_evaluator.cc
 )

file(GLOB ats_eos_inc_files "*.hh")

set(ats_eos_link_libs
  ${Teuchos_LIBRARIES}
  ${Epetra_LIBRARIES}
  error_handling
  atk
  mesh
  data_structures
  whetstone
  solvers
  state
  )


add_amanzi_library(ats_eos
                   SOURCE ${ats_eos_src_files}
                   HEADERS ${ats_eos_inc_files}
		   LINK_LIBS ${ats_eos_link_libs})

