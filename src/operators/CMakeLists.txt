# -*- mode: cmake -*-

#
#  ATS Operators
#
include_directories(${ATS_SOURCE_DIR}/src/operators/advection)
include_directories(${ATS_SOURCE_DIR}/src/operators/upwinding)
include_directories(${ATS_SOURCE_DIR}/src/operators/deformation)

set(ats_operators_src_files
  advection/advection.cc
  advection/advection_donor_upwind.cc
  advection/advection_factory.cc
  upwinding/upwind_cell_centered.cc
  upwinding/upwind_arithmetic_mean.cc
  upwinding/UpwindFluxFactory.cc
  upwinding/upwind_flux_split_denominator.cc  
  upwinding/upwind_flux_fo_cont.cc 
  upwinding/upwind_flux_harmonic_mean.cc
  upwinding/upwind_elevation_stabilized.cc
  upwinding/upwind_total_flux.cc
  upwinding/upwind_potential_difference.cc
  upwinding/upwind_gravity_flux.cc
#  deformation/MatrixVolumetricDeformation.cc
#  deformation/Matrix_PreconditionerDelegate.cc
  )

set(ats_operators_inc_files
  advection/advection.hh
  advection/advection_donor_upwind.hh
  advection/advection_factory.hh
  upwinding/upwinding.hh
  upwinding/UpwindFluxFactory.hh
  upwinding/upwind_arithmetic_mean.hh
  upwinding/upwind_cell_centered.hh
  upwinding/upwind_flux_fo_cont.hh
  upwinding/upwind_flux_harmonic_mean.hh
  upwinding/upwind_flux_split_denominator.hh
  upwinding/upwind_gravity_flux.hh
  upwinding/upwind_potential_difference.hh
  upwinding/upwind_elevation_stabilized.hh
  upwinding/upwind_total_flux.hh
#  deformation/MatrixVolumetricDeformation.hh
#  deformation/Matrix_PreconditionerDelegate.hh
  )


set(ats_operators_link_libs
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


add_amanzi_library(ats_operators
                   SOURCE ${ats_operators_src_files}
                   HEADERS ${ats_operators_inc_files}
		   LINK_LIBS ${ats_operators_link_libs})

