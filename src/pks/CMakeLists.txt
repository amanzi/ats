# -*- mode: cmake -*-

#
#  ATS
#    PK class
#

set(ats_pks_src_files
  pk_helpers.cc
  pk_bdf_default.cc
  pk_physical_default.cc
  pk_physical_bdf_default.cc
  pk_explicit_default.cc
  bc_factory.cc
  )

file(GLOB ats_pks_inc_files "*.hh")

set(ats_pks_link_libs
  ${Teuchos_LIBRARIES}
  ${Epetra_LIBRARIES}
  error_handling
  atk
  mesh
  data_structures
  whetstone
  solvers
  state
  time_integration
  pks
  )


add_amanzi_library(ats_pks
                   SOURCE ${ats_pks_src_files}
                   HEADERS ${ats_pks_inc_files}
		   LINK_LIBS ${ats_pks_link_libs})


add_subdirectory(energy)
add_subdirectory(flow)
add_subdirectory(transport)
add_subdirectory(deform)
add_subdirectory(surface_balance)
add_subdirectory(biogeochemistry)
add_subdirectory(mpc)
