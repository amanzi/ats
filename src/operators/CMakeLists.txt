# -*- mode: cmake -*-

#
#  ATS Operators
#
# include_directories(${ATS_SOURCE_DIR}/src/operators/advection)
# include_directories(${ATS_SOURCE_DIR}/src/operators/deformation)

set(ats_operators_src_files
  # advection/advection.cc
  # advection/advection_donor_upwind.cc
  # advection/advection_factory.cc
#  deformation/MatrixVolumetricDeformation.cc
#  deformation/Matrix_PreconditionerDelegate.cc
  )

set(ats_operators_inc_files
  # advection/advection.hh
  # advection/advection_donor_upwind.hh
  # advection/advection_factory.hh
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

