# -*- mode: cmake -*-

#
#  ATS
#    Transport Amanzi
#

#================================================
add_subdirectory(constitutive_relations)

# ATS include directories
include_directories(${ATS_SOURCE_DIR}/src/pks)
include_directories(${ATS_SOURCE_DIR}/src/pks/transport/constitutive_relations/sources)
# include_directories(${ATS_SOURCE_DIR}/src/pks/transport/constitutive_relations/sediment_transport)

# Amanzi include directories
include_directories(${TRANSPORT_SOURCE_DIR})

set(ats_transport_src_files
  #transport_ats_dispersion.cc
  transport_ats_ti.cc
  #transport_ats_henrylaw.cc
  #transport_ats_vandv.cc
  transport_ats_pk.cc
  # sediment_transport_pk.cc
 )


set(ats_transport_inc_files
  transport_ats.hh
  # sediment_transport_pk.hh
  )


set(ats_transport_link_libs
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
  chemistry_pk
  transport
  ats_operators
  ats_eos
  ats_pks
  )


add_amanzi_library(ats_transport
                   SOURCE ${ats_transport_src_files}
                   HEADERS ${ats_transport_inc_files}
		   LINK_LIBS ${ats_transport_link_libs})

#================================================
# register evaluators/factories/pks

# transport_amanzi
register_evaluator_with_factory(
  HEADERFILE transport_ats_reg.hh
  LISTNAME   ATS_TRANSPORT_REG
  )

# register_evaluator_with_factory(
#   HEADERFILE sediment_transport_reg.hh
#   LISTNAME   ATS_TRANSPORT_REG
#   )
  
generate_evaluators_registration_header(
  HEADERFILE ats_transport_registration.hh
  LISTNAME   ATS_TRANSPORT_REG
  INSTALL    True
  )

