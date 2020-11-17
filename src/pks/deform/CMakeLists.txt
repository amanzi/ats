# -*- mode: cmake -*-

#
#  ATS
#    Deformation PK class
#
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/constitutive_relations/porosity)
include_directories(${ATS_SOURCE_DIR}/src/operators/deformation)
include_directories(${ATS_SOURCE_DIR}/src/pks)

set(ats_deform_src_files
  constitutive_relations/porosity/porosity_evaluator.cc
  volumetric_deformation.cc
  )

set(ats_deform_inc_files
  constitutive_relations/porosity/porosity_evaluator.hh
  volumetric_deformation.hh
  )


set(ats_deform_link_libs
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
  ats_pks
  )


add_amanzi_library(ats_deform
                   SOURCE ${ats_deform_src_files}
                   HEADERS ${ats_deform_inc_files}
		   LINK_LIBS ${ats_deform_link_libs})


# volumetric_deformation/
register_evaluator_with_factory(
  HEADERFILE volumetric_deformation_reg.hh
  LISTNAME  ATS_DEFORMATION_REG)

register_evaluator_with_factory(
  HEADERFILE constitutive_relations/porosity/porosity_evaluator_reg.hh
  LISTNAME  ATS_DEFORMATION_REG)
generate_evaluators_registration_header(
  HEADERFILE ats_deformation_registration.hh
  LISTNAME   ATS_DEFORMATION_REG
  INSTALL    True
  )


  
