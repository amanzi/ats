# -*- mode: cmake -*-
#
#  ATS
#    Constitutive relations for transport
#


# collect all sources
list(APPEND subdirs sources)
#list(APPEND subdirs sources sediment_transport)

set(ats_transport_relations_src_files "")
set(ats_transport_relations_inc_files "")

foreach(lcv IN LISTS subdirs)
  include_directories(${ATS_SOURCE_DIR}/src/pks/transport/constitutive_relations/${lcv})
  
  file(GLOB subdir_sources "./${lcv}/*.cc")
  set(ats_transport_relations_src_files ${ats_transport_relations_src_files} ${subdir_sources})

  file(GLOB subdir_incs "./${lcv}/*.hh")
  set(ats_transport_relations_inc_files ${ats_transport_relations_inc_files} ${subdir_incs})

  file(GLOB registrations "./${lcv}/*_reg.hh" )
  foreach(reg_lcv IN LISTS registrations)
    register_abs_evaluator_with_factory(HEADERFILE ${reg_lcv} LISTNAME ATS_TRANSPORT_RELATIONS_REG)
  endforeach(reg_lcv)

endforeach(lcv)

set(ats_transport_relations_link_libs
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
add_amanzi_library(ats_transport_relations
                   SOURCE ${ats_transport_relations_src_files}
                   HEADERS ${ats_transport_relations_inc_files}
		   LINK_LIBS ${ats_transport_relations_link_libs})

generate_evaluators_registration_header(
    HEADERFILE ats_transport_relations_registration.hh
    LISTNAME   ATS_TRANSPORT_RELATIONS_REG
    INSTALL    True
    )
                 
