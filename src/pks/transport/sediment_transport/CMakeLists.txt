# -*- mode: cmake -*-

#
#  
#   Sediment transport process kernel
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
##include(PrintVariable)
##include(TestManager)

#
# Define a project name
# After this command the following varaibles are defined
#   TRANSPORT_SOURCE_DIR
#   TRANSPORT_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
##project(TRANSPORT)


# External (TPL) include directories
##include_directories(${Epetra_INCLUDE_DIRS})
##include_directories(${ASCEMIO_INCLUDE_DIR})
# Amanzi include directories
##include_directories(${ATK_SOURCE_DIR})
##include_directories(${CHEM_INCLUDES_DIR})
##include_directories(${CHEMPK_SOURCE_DIR})
##include_directories(${DATA_STRUCTURES_SOURCE_DIR})
##include_directories(${FUNCS_SOURCE_DIR})
##include_directories(${GEOMETRY_SOURCE_DIR})
##include_directories(${MESH_SOURCE_DIR})
##include_directories(${MESH_FACTORY_SOURCE_DIR})
##include_directories(${MFUNCS_SOURCE_DIR})
##include_directories(${OUTPUT_SOURCE_DIR})
##include_directories(${OPERATORS_SOURCE_DIR})
##include_directories(${PKS_SOURCE_DIR})
##include_directories(${SOLVERS_SOURCE_DIR})
##include_directories(${STATE_SOURCE_DIR})
##include_directories(${TIME_INTEGRATION_SOURCE_DIR})
##include_directories(${WHETSTONE_SOURCE_DIR})
include_directories(${FUNCTIONS_SOURCE_DIR})
include_directories(${TRANSPORT_SOURCE_DIR})
include_directories(${ATS_SOURCE_DIR}/src/pks)
include_directories(${ATS_SOURCE_DIR}/src/pks/transport/)


#include_directories(${Amanzi_TPL_MSTK_INCLUDE_DIRS})

#
# Transport registrations
#



#
# Library: transport
#
set(sed_transport_inc_files
  sediment_transport_pk.hh
  erosion_evaluator.hh
  settlement_evaluator.hh
  trapping_evaluator.hh
  organic_matter_evaluator.hh
  )

set(sed_transport_src_files
  sediment_transport_pk.cc
  Transport_Diffusion.cc
  erosion_evaluator.cc
  settlement_evaluator.cc
  trapping_evaluator.cc
  organic_matter_evaluator.cc)

set(sed_transport_link_libs
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
    )

add_amanzi_library(ats_sed_transport
                   SOURCE ${sed_transport_src_files}
                   HEADERS ${sed_transport_inc_files}
		   LINK_LIBS ${sed_transport_link_libs})

register_evaluator_with_factory(
  HEADERFILE sediment_transport_reg.hh
  LISTNAME   SED_TRANSPORT_REG
  )


#sediment diposition and erosion rates
register_evaluator_with_factory(
  HEADERFILE erosion_evaluator_reg.hh
  LISTNAME SED_TRANSPORT_REG
  )
register_evaluator_with_factory(
  HEADERFILE settlement_evaluator_reg.hh
  LISTNAME SED_TRANSPORT_REG
  )
register_evaluator_with_factory(
  HEADERFILE trapping_evaluator_reg.hh
  LISTNAME SED_TRANSPORT_REG
  )
register_evaluator_with_factory(
  HEADERFILE organic_matter_evaluator_reg.hh
  LISTNAME SED_TRANSPORT_REG
  )

generate_evaluators_registration_header( 
  HEADERFILE ats_sediment_transport_registration.hh   
  LISTNAME SED_TRANSPORT_REG  
  INSTALL  True 
  )
