# -*- mode: cmake -*-

#
#  ATS
#    Transport Amanzi
#

#================================================
# register evaluators/factories/pks

include_directories(${CHEMPK_SOURCE_DIR})
include_directories(${AMANZI_SOURCE_DIR}/src/common/alquimia)
include_directories(${FUNCTIONS_SOURCE_DIR})
include_directories(${TRANSPORT_SOURCE_DIR})
include_directories(${ATS_SOURCE_DIR}/src/pks)

set(ats_transport_src_files
  transport_ats_dispersion.cc
  transport_ats_ti.cc
  transport_ats_henrylaw.cc
  transport_ats_vandv.cc
  transport_ats_initialize.cc
  transport_ats_pk.cc
  # sediment_transport/sediment_transport_pk.cc
  # sediment_transport/erosion_evaluator.cc
  # sediment_transport/settlement_evaluator.cc
  # sediment_transport/trapping_evaluator.cc
  # sediment_transport/organic_matter_evaluator.cc  
 )


set(ats_transport_inc_files
  transport_ats.hh
  # sediment_transport/sediment_transport_pk.hh
  # sediment_transport/erosion_evaluator.hh
  # sediment_transport/settlement_evaluator.hh
  # sediment_transport/trapping_evaluator.hh
  # sediment_transport/organic_matter_evaluator.hh
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

if (APPLE)
  set_target_properties(ats_transport PROPERTIES LINK_FLAGS "-ld_classic")
endif()

#================================================
# register evaluators/factories/pks

# transport_amanzi
register_evaluator_with_factory(
  HEADERFILE transport_ats_reg.hh
  LISTNAME   ATS_TRANSPORT_REG
  )

# register_evaluator_with_factory(
#   HEADERFILE sediment_transport/sediment_transport_reg.hh
#   LISTNAME   ATS_TRANSPORT_REG
#   )
  
# # sediment diposition and erosion rates
# register_evaluator_with_factory(
#   HEADERFILE sediment_transport/erosion_evaluator_reg.hh
#   LISTNAME ATS_TRANSPORT_REG
#   )
# register_evaluator_with_factory(
#   HEADERFILE sediment_transport/settlement_evaluator_reg.hh
#   LISTNAME ATS_TRANSPORT_REG
#   )
# register_evaluator_with_factory(
#   HEADERFILE sediment_transport/trapping_evaluator_reg.hh
#   LISTNAME ATS_TRANSPORT_REG
#   )
# register_evaluator_with_factory(
#   HEADERFILE sediment_transport/organic_matter_evaluator_reg.hh
#   LISTNAME ATS_TRANSPORT_REG
#   )

  
generate_evaluators_registration_header(
  HEADERFILE ats_transport_registration.hh
  LISTNAME   ATS_TRANSPORT_REG
  INSTALL    True
  )

