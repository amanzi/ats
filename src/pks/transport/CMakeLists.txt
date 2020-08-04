# -*- mode: cmake -*-

#
#  ATS
#    Transport Amanzi
#


# ATS include directories
add_subdirectory(sediment_transport)

#================================================
# register evaluators/factories/pks

include_directories(${CHEMPK_SOURCE_DIR})
include_directories(${AMANZI_SOURCE_DIR}/src/common/alquimia)
include_directories(${FUNCTIONS_SOURCE_DIR})


set(ats_transport_src_files
  MDM_BurnettFrind.cc
  MDM_Bear.cc
  MDMPartition.cc
  MDM_LichtnerKelkarRobinson.cc
  MDM_Isotropic.cc
  MDMFactory.cc
  MultiscaleTransportPorosity_DPM.cc
  Transport_Dispersion.cc
  Transport_TI.cc
  MultiscaleTransportPorosity_GDPM.cc
  Transport_HenryLaw.cc
  Transport_VandV.cc
  MultiscaleTransportPorosityFactory.cc
  TransportBoundaryFunction_Alquimia.cc
  Transport_Initialize.cc
  MultiscaleTransportPorosityPartition.cc
  TransportSourceFunction_Alquimia.cc
  Transport_PK.cc
 )


set(ats_transport_inc_files
  TransportDefs.hh
  MDM.hh
  MDM_Isotropic.hh
  MDMFactory.hh
  MDM_LichtnerKelkarRobinson.hh
  MDMPartition.hh
  MDM_Bear.hh
  MDM_BurnettFrind.hh
  DiffusionPhase.hh
  MaterialProperties.hh
  MultiscaleTransportPorosityPartition.hh
  TransportDomainFunction.hh  
  MultiscaleTransportPorosity_DPM.hh
  TransportSourceFunction_Alquimia.hh
  MultiscaleTransportPorosity_GDPM.hh
  MultiscaleTransportPorosity.hh
  MultiscaleTransportPorosityFactory.hh
  TransportBoundaryFunction_Alquimia.hh
  Transport_PK_ATS.hh
  )


set(ats_transport_link_libs
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
  chemistry_pk
  ats_operators
  ats_eos
  ats_pks
  )


add_amanzi_library(ats_transport
                   SOURCE ${ats_transport_src_files}
                   HEADERS ${ats_transport_inc_files}
		   LINK_LIBS ${ats_transport_link_libs})

#================================================
# register evaluators/factories/pks

# transport_amanzi
register_evaluator_with_factory(
  HEADERFILE Transport_PK_reg.hh
  LISTNAME   ATS_TRANSPORT_REG
  )

# -- mechamical dispersion models registration
register_evaluator_with_factory(
  HEADERFILE MDMFactory_reg.hh
  LISTNAME   ATS_TRANSPORT_REG
  )
register_evaluator_with_factory(
  HEADERFILE MDM_Isotropic_reg.hh
  LISTNAME   ATS_TRANSPORT_REG
  )
register_evaluator_with_factory(
  HEADERFILE MDM_Bear_reg.hh
  LISTNAME   ATS_TRANSPORT_REG
  )
register_evaluator_with_factory(
  HEADERFILE MDM_BurnettFrind_reg.hh
  LISTNAME   ATS_TRANSPORT_REG
  )
register_evaluator_with_factory(
  HEADERFILE MDM_LichtnerKelkarRobinson_reg.hh
  LISTNAME   ATS_TRANSPORT_REG
  )

# -- multiscale models registration
register_evaluator_with_factory(
  HEADERFILE MultiscaleTransportPorosityFactory_reg.hh
  LISTNAME   ATS_TRANSPORT_REG
  )
register_evaluator_with_factory(
  HEADERFILE MultiscaleTransportPorosity_DPM_reg.hh
  LISTNAME   ATS_TRANSPORT_REG
  )
register_evaluator_with_factory(
  HEADERFILE MultiscaleTransportPorosity_GDPM_reg.hh
  LISTNAME   ATS_TRANSPORT_REG
  )

generate_evaluators_registration_header(
  HEADERFILE ats_transport_registration.hh
  LISTNAME   ATS_TRANSPORT_REG
  INSTALL    True
  )

