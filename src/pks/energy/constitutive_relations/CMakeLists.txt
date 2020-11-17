# -*- mode: cmake -*-

#
#  ATS
#    Constitutive relations for flow
#

# collect all registrations and register
file(GLOB_RECURSE registrations "./*_reg.hh" )
foreach(reg_lcv IN LISTS registrations)
  register_abs_evaluator_with_factory(HEADERFILE ${reg_lcv} LISTNAME ATS_ENERGY_RELATIONS_REG)
endforeach(reg_lcv)

generate_evaluators_registration_header(
  HEADERFILE ats_energy_relations_registration.hh
  LISTNAME   ATS_ENERGY_RELATIONS_REG
  INSTALL    True
  )

# collect all sources
list(APPEND subdirs energy enthalpy internal_energy source_terms thermal_conductivity)
set(ats_energy_relations_src_files "")
set(ats_energy_relations_inc_files "")

foreach(lcv IN LISTS subdirs)
  include_directories(${ATS_SOURCE_DIR}/src/pks/energy/constitutive_relations/${lcv})
  
  file(GLOB subdir_sources "./${lcv}/*.cc")
  set(ats_energy_relations_src_files ${ats_energy_relations_src_files} ${subdir_sources})

  file(GLOB subdir_incs "./${lcv}/*.hh")
  set(ats_energy_relations_inc_files ${ats_energy_relations_inc_files} ${subdir_incs})
endforeach(lcv)

set(ats_energy_relations_link_libs
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
add_amanzi_library(ats_energy_relations
                   SOURCE ${ats_energy_relations_src_files}
                   HEADERS ${ats_energy_relations_inc_files}
		   LINK_LIBS ${ats_energy_relations_link_libs})


