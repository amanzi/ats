# -*- mode: cmake -*-

#
#  ATS
#    Constitutive relations for lake thermo
#

# collect all registrations and register
file(GLOB_RECURSE registrations "./*_reg.hh" )
foreach(reg_lcv IN LISTS registrations)
  register_abs_evaluator_with_factory(HEADERFILE ${reg_lcv} LISTNAME ATS_LAKE_THERMO_RELATIONS_REG)
endforeach(reg_lcv)

generate_evaluators_registration_header(
  HEADERFILE ats_lake_thermo_relations_registration.hh
  LISTNAME   ATS_LAKE_THERMO_RELATIONS_REG
  INSTALL    True
  )

# collect all sources
list(APPEND subdirs density energy enthalpy internal_energy thermal_conductivity heat_capacity heat_flux_bc evaporation_rate surface_temperature)
set(ats_lake_thermo_relations_src_files "")
set(ats_lake_thermo_relations_inc_files "")

foreach(lcv IN LISTS subdirs)
  include_directories(${ATS_SOURCE_DIR}/src/pks/lake_thermo/constitutive_relations/${lcv})
  
  file(GLOB subdir_sources "./${lcv}/*.cc")
  set(ats_lake_thermo_relations_src_files ${ats_lake_thermo_relations_src_files} ${subdir_sources})

  file(GLOB subdir_incs "./${lcv}/*.hh")
  set(ats_lake_thermo_relations_inc_files ${ats_lake_thermo_relations_inc_files} ${subdir_incs})
endforeach(lcv)

set(ats_lake_thermo_relations_link_libs
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

# make the library
add_amanzi_library(ats_lake_thermo_relations
                   SOURCE ${ats_lake_thermo_relations_src_files}
                   HEADERS ${ats_lake_thermo_relations_inc_files}
		   LINK_LIBS ${ats_lake_thermo_relations_link_libs})


