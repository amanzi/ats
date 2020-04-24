# -*- mode: cmake -*-


# generic PKs
include_directories(${ATS_SOURCE_DIR}/pks)


set(ats_basic_inc_files
  ConservationEquation.hh
  ConservationODE.hh
  MixedFormParabolicPDE.hh
  )  

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
