# -*- mode: cmake -*-

#
#  ATS
#    Lake thermo PK class
#

# additional common components
add_subdirectory(constitutive_relations)

# ATS include directories
include_directories(${ATS_SOURCE_DIR}/src/pks)
include_directories(${ATS_SOURCE_DIR}/src/operators/advection)
include_directories(${ATS_SOURCE_DIR}/src/operators/upwinding)
include_directories(${ATS_SOURCE_DIR}/src/pks/lake_thermo/constitutive_relations/density)
include_directories(${ATS_SOURCE_DIR}/src/pks/lake_thermo/constitutive_relations/energy)
include_directories(${ATS_SOURCE_DIR}/src/pks/lake_thermo/constitutive_relations/enthalpy)
include_directories(${ATS_SOURCE_DIR}/src/pks/lake_thermo/constitutive_relations/thermal_conductivity)
include_directories(${ATS_SOURCE_DIR}/src/pks/lake_thermo/constitutive_relations/heat_capacity)
include_directories(${ATS_SOURCE_DIR}/src/pks/lake_thermo/constitutive_relations/heat_flux_bc)
include_directories(${ATS_SOURCE_DIR}/src/pks/lake_thermo/constitutive_relations/evaporation_rate)
include_directories(${ATS_SOURCE_DIR}/src/pks/lake_thermo/constitutive_relations/surface_temperature)
include_directories(${ATS_SOURCE_DIR}/src/constitutive_relations/eos)


set(ats_lake_src_files
  lake_thermo_pk.cc
  lake_thermo_physics.cc
  lake_thermo_ti.cc
  )

set(ats_lake_inc_files
  lake_thermo_pk.hh
  lake_thermo_bc_factory.hh
  )


set(ats_lake_link_libs
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
  ats_lake_thermo_relations
  )


add_amanzi_library(ats_lake
                   SOURCE ${ats_lake_src_files}
                   HEADERS ${ats_lake_inc_files}
		   LINK_LIBS ${ats_lake_link_libs})



#================================================
# register evaluators/factories/pks
register_evaluator_with_factory(
  HEADERFILE lake_thermo_pk_reg.hh
  LISTNAME   ATS_LAKE_THERMO_PKS_REG
  )

generate_evaluators_registration_header(
  HEADERFILE ats_lake_pks_registration.hh
  LISTNAME   ATS_LAKE_THERMO_PKS_REG
  INSTALL    True
  )

