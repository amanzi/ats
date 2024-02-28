/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

#pragma once

#include "EvaluatorModelCV.hh"
#include "EvaluatorModelCVByMaterial.hh"
#include "Factory.hh"

#include "registration_macro.hh"
#include "capillary_pressure_liquid_atm_model.hh"

#include "wrm_van_genuchten.hh"
#include "wrm_model.hh"
#include "relative_permeability_model.hh"
#include "relative_hydraulic_conductivity_evaluator.hh"

namespace Amanzi {

REGISTER_MODEL(Flow::Relations::CapillaryPressureLiquidAtmModel);

REGISTER_MODEL(Flow::Relations::WRMVanGenuchtenModel);
REGISTER_BY_MATERIAL(Flow::Relations::WRMVanGenuchtenModel);

REGISTER_MODEL(Flow::Relations::RelativePermeabilityVanGenuchtenModel);
REGISTER_BY_MATERIAL(Flow::Relations::RelativePermeabilityVanGenuchtenModel);

REGISTER(Flow::Relations::RelativeHydraulicConductivityEvaluator);

} // namespace Amanzi
