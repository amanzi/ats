# -*- mode: cmake -*-


# generic PKs
include_directories(${ATS_SOURCE_DIR}/pks)


set(ats_basic_inc_files
  ConservationEquation.hh
  ConservationODE.hh
  )  


#================================================
# register evaluators/factories/pks

register_evaluator_with_factory(
  HEADERFILE ConservationODE_reg.hh
  LISTNAME ATS_REG
)
