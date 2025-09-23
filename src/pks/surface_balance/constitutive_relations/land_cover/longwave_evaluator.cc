/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "Key.hh"
#include "seb_physics_defs.hh"
#include "seb_physics_funcs.hh"
#include "longwave_evaluator.hh"


namespace Amanzi {
namespace ATS_Physics {
namespace SurfaceBalance {
namespace Relations {

LongwaveEvaluator::LongwaveEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  auto domain = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  air_temp_key_ = Keys::readKey(plist, domain, "air temperature", "air_temperature");
  dependencies_.insert(KeyTag{ air_temp_key_, tag });
  vp_air_key_ = Keys::readKey(plist, domain, "vapor pressure air", "vapor_pressure_air");
  dependencies_.insert(KeyTag{ vp_air_key_, tag });

  scale_ = plist.get<double>("scaling factor [-]", 1.0);
}

// Required methods from EvaluatorSecondaryMonotypeCV
void
LongwaveEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  const auto& air_temp = *S.Get<CompositeVector>(air_temp_key_, tag).ViewComponent("cell", false);
  const auto& vp_air = *S.Get<CompositeVector>(vp_air_key_, tag).ViewComponent("cell", false);
  auto& res = *result[0]->ViewComponent("cell", false);

  for (int c = 0; c != res.MyLength(); ++c) {
    res[0][c] = scale_ * Relations::IncomingLongwaveRadiation(air_temp[0][c], vp_air[0][c]);
  }
}

} // namespace Relations
} // namespace SurfaceBalance
} // namespace ATS_Physics
} // namespace Amanzi
