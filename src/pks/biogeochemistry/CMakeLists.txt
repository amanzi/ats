# -*- mode: cmake -*-
include_directories(${ATS_SOURCE_DIR}/src/pks)

# ATS Surface balance PKs describe Evaporation, energy fluxes from
#  long/showtwave radiation, precip, etc etc etc


if (ENABLE_FATES)
  add_subdirectory(fates)
endif()

set(ats_bgc_src_files
  bgc_simple/utils.cc
  bgc_simple/PFT.cc
  bgc_simple/SoilCarbonParameters.cc
  bgc_simple/vegetation.cc
  bgc_simple/bgc_simple_funcs.cc
  bgc_simple/bgc_simple.cc
  carbon/simple/CarbonSimple.cc
  constitutive_models/carbon/bioturbation_evaluator.cc  
  )

set(ats_bgc_inc_files
  bgc_simple/utils.hh
  bgc_simple/PFT.hh
  bgc_simple/SoilCarbonParameters.hh
  bgc_simple/vegetation.hh
  bgc_simple/bgc_simple_funcs.hh
  bgc_simple/bgc_simple.hh
  carbon/simple/CarbonSimple.hh
  constitutive_models/carbon/bioturbation_evaluator.hh
  )


set(ats_bgc_link_libs
  ${Teuchos_LIBRARIES}
  ${Epetra_LIBRARIES}
  error_handling
  atk
  mesh
  data_structures
  state
  pks
  ats_pks
  )


add_amanzi_library(ats_bgc
                   SOURCE ${ats_bgc_src_files}
                   HEADERS ${ats_bgc_inc_files}
		   LINK_LIBS ${ats_bgc_link_libs})

#================================================
# register evaluators/factories/pks


register_evaluator_with_factory(
  HEADERFILE bgc_simple/bgc_simple_reg.hh
  LISTNAME ATS_BGC_REG
)

if (ENABLE_FATES)
  register_evaluator_with_factory(
    HEADERFILE fates/fates_reg.hh
    LISTNAME ATS_BGC_REG
  )
endif()

register_evaluator_with_factory(
  HEADERFILE carbon/simple/CarbonSimple_reg.hh
  LISTNAME ATS_BGC_REG
)

register_evaluator_with_factory(
  HEADERFILE constitutive_models/carbon/bioturbation_evaluator_reg.hh
  LISTNAME ATS_BGC_REG
)

generate_evaluators_registration_header(
  HEADERFILE ats_bgc_registration.hh
  LISTNAME   ATS_BGC_REG
  INSTALL    True
)


