# -*- mode: cmake -*-

#
#  ATS
#    Energy-water content globalizations.
#
#

include_directories(${ATS_SOURCE_DIR}/src/constitutive_relations/eos)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/wrm)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/elevation)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/porosity)
include_directories(${ATS_SOURCE_DIR}/src/pks/energy/constitutive_relations/internal_energy)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/overland_conductivity)

set(ats_mpc_relations_src_files
  ewc_model_base.cc
  liquid_ice_model.cc
  permafrost_model.cc
  surface_ice_model.cc
 )

set(ats_mpc_relations_inc_files
  ewc_model.hh
  ewc_model_base.hh
  liquid_ice_model.hh
  permafrost_model.hh
  surface_ice_model.hh
  thermal_richards_model.hh
  )


set(ats_mpc_relations_link_libs
  ${Teuchos_LIBRARIES}
  ${Epetra_LIBRARIES}
  error_handling
  atk
  mesh
  data_structures
  whetstone
  solvers
  state
  ats_eos
  ats_flow_relations
  ats_energy_relations
  )


add_amanzi_library(ats_mpc_relations
                   SOURCE ${ats_mpc_relations_src_files}
                   HEADERS ${ats_mpc_relations_inc_files}
		   LINK_LIBS ${ats_mpc_relations_link_libs})





