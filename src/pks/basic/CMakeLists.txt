# -*- mode: cmake -*-

# generic PKs
include_directories(
  ${ATS_SOURCE_DIR}/pks
  ${STATE_SOURCE_DIR}
  ${STATE_SOURCE_DIR}/evaluator
  ${STATE_SOURCE_DIR}/data
  ${STATE_SOURCE_DIR}/io
  )


set(ats_basic_inc_files
  PK_MixinConservationEquation.hh
  ConservationODE.hh
  MixedFormParabolicPDE.hh
  EvaluatorUpwindPotentialDifference.hh
  )

set(ats_basic_src_files
  EvaluatorUpwindPotentialDifference.cc
  )

set(ats_basic_link_libs
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
  )


add_amanzi_library(ats_basic
                   SOURCE ${ats_basic_src_files}
                   HEADERS ${ats_basic_inc_files}
		   LINK_LIBS ${ats_basic_link_libs})
                 

#================================================
# register evaluators/factories/pks
register_evaluator_with_factory(
  HEADERFILE ConservationODE_reg.hh
  LISTNAME ATS_REG
)

register_evaluator_with_factory(
  HEADERFILE MixedFormParabolicPDE_reg.hh
  LISTNAME ATS_REG
)

register_evaluator_with_factory(
  HEADERFILE EvaluatorUpwindPotentialDifference_reg.hh
  LISTNAME ATS_REG
)

