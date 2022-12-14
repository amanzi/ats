/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "state_evaluators_registration.hh"

#include "ats_relations_registration.hh"
#include "ats_transport_registration.hh"
#include "ats_energy_pks_registration.hh"
#include "ats_energy_relations_registration.hh"
#include "ats_flow_pks_registration.hh"
#include "ats_flow_relations_registration.hh"
#include "ats_deformation_registration.hh"
#include "ats_bgc_registration.hh"
#include "ats_surface_balance_registration.hh"
#include "ats_mpc_registration.hh"
//#include "ats_sediment_transport_registration.hh"
#include "mdm_transport_registration.hh"
#include "multiscale_transport_registration.hh"
#ifdef ALQUIMIA_ENABLED
#  include "pks_chemistry_registration.hh"
#endif
