# -*- mode: cmake -*-

project(ATS)
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${DBC_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})
include_directories(${OUTPUT_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${FUNCS_SOURCE_DIR})
include_directories(${MFUNCS_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${DATA_STRUCTURES_SOURCE_DIR})
include_directories(${STATE_SOURCE_DIR})
include_directories(${WHETSTONE_SOURCE_DIR})
include_directories(${OPERATORS_SOURCE_DIR})
include_directories(${SOLVERS_SOURCE_DIR})
include_directories(${TIME_INTEGRATION_SOURCE_DIR})
include_directories(${PKS_SOURCE_DIR})

# # operators -- layer between discretization and PK
# add_subdirectory(operators)

# add_subdirectory(constitutive_relations)

# Process Kernels
add_subdirectory(pks)


generate_evaluators_registration_header(
  HEADERFILE executables/ats_registration.hh
  LISTNAME   ATS_REG
  INSTALL    True
  )


# main
add_subdirectory(executables)
