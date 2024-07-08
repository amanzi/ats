#include "EvaluatorModelCV.hh"
#include "EvaluatorModelCVByMaterial.hh"
#include "Factory.hh"

#include "registration_macro.hh"
#include "richards_water_content_model.hh"
#include "overland_pressure_water_content_model.hh"
#include "ponded_depth_model.hh"

namespace Amanzi {

REGISTER_MODEL(Flow::Relations::RichardsWaterContentModel);
REGISTER_MODEL(Flow::Relations::OverlandPressureWaterContentModel);
REGISTER_MODEL(Flow::Relations::PondedDepthModel);

} // namespace Amanzi
#include "EvaluatorModelCV.hh"
#include "EvaluatorModelCVByMaterial.hh"
#include "Factory.hh"

#include "registration_macro.hh"
#include "compressible_porosity_linear_model.hh"

namespace Amanzi {

REGISTER_MODEL(Flow::Relations::CompressiblePorosityLinearModel);
REGISTER_BY_MATERIAL(Flow::Relations::CompressiblePorosityLinearModel);

} // namespace Amanzi
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
#include "EvaluatorModelCV.hh"
#include "EvaluatorModelCVByMaterial.hh"
#include "Factory.hh"

#include "registration_macro.hh"
#include "soil_resistance_sakagucki_zeng_evaluator.hh"

namespace Amanzi {

template <>
REGISTER(Flow::Relations::SoilResistanceSakaguckiZengEvaluator);
template <>
REGISTER(Flow::Relations::SoilResistanceSakaguckiZengEvaluatorByMaterial);

} // namespace Amanzi
#include "Factory.hh"

#include "registration_macro.hh"

#include "overland_conductivity_evaluator.hh"

namespace Amanzi {

REGISTER(Flow::Relations::OverlandConductivityEvaluator);


} // namespace Amanzi
