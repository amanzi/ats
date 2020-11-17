# -*- mode: cmake -*-

#
#  ATS
#    Flow PK class
#

# additional common components
add_subdirectory(constitutive_relations)

# ATS include directories
include_directories(${ATS_SOURCE_DIR}/src/pks)
include_directories(${ATS_SOURCE_DIR}/src/operators/advection)
include_directories(${ATS_SOURCE_DIR}/src/operators/upwinding)
include_directories(${ATS_SOURCE_DIR}/src/pks/energy/constitutive_relations/enthalpy)
include_directories(${ATS_SOURCE_DIR}/src/pks/energy/constitutive_relations/energy)
include_directories(${ATS_SOURCE_DIR}/src/pks/energy/constitutive_relations/internal_energy)
include_directories(${ATS_SOURCE_DIR}/src/pks/energy/constitutive_relations/thermal_conductivity)
include_directories(${ATS_SOURCE_DIR}/src/pks/energy/constitutive_relations/source_terms)
include_directories(${ATS_SOURCE_DIR}/src/constitutive_relations/eos)


set(ats_energy_src_files
  energy_base_pk.cc
  energy_base_physics.cc
  energy_base_ti.cc
  energy_surface_ice.cc
  energy_two_phase.cc
  energy_three_phase.cc
  energy_interfrost.cc
  )

set(ats_energy_inc_files
  energy_bc_factory.hh
  energy_base.hh
  energy_surface_ice.hh
  energy_two_phase.hh
  energy_three_phase.hh
  energy_interfrost.hh
  )


set(ats_energy_link_libs
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
  ats_eos
  ats_pks
  ats_energy_relations
  )


add_amanzi_library(ats_energy
                   SOURCE ${ats_energy_src_files}
                   HEADERS ${ats_energy_inc_files}
		   LINK_LIBS ${ats_energy_link_libs})



#================================================
# register evaluators/factories/pks
register_evaluator_with_factory(
  HEADERFILE energy_surface_ice_reg.hh
  LISTNAME   ATS_ENERGY_PKS_REG
  )

register_evaluator_with_factory(
  HEADERFILE energy_three_phase_reg.hh
  LISTNAME   ATS_ENERGY_PKS_REG
  )

register_evaluator_with_factory(
  HEADERFILE energy_two_phase_reg.hh
  LISTNAME   ATS_ENERGY_PKS_REG
  )

generate_evaluators_registration_header(
  HEADERFILE ats_energy_pks_registration.hh
  LISTNAME   ATS_ENERGY_PKS_REG
  INSTALL    True
  )

