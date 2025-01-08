/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ahmad Jan (jana@ornl.gov)
*/

/*
  The elevation evaluator gets the subsurface temperature and computes moisture content
  over time.

*/

#include "carbon_decomposition_rate_evaluator.hh"

namespace Amanzi {
namespace Relations {


CarbonDecomposeRateEvaluator::CarbonDecomposeRateEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  Tag tag = my_keys_.front().second;
  Key domain = Keys::getDomain(my_keys_.front().first); //column domain

  temp_key_ = Keys::readKey(plist, domain, "temperature", "temperature");
  dependencies_.insert(KeyTag{ temp_key_, tag });

  pres_key_ = Keys::readKey(plist, domain, "pressure", "pressure");
  dependencies_.insert(KeyTag{ pres_key_, tag });

  por_key_ = Keys::readKey(plist, domain, "porosity", "porosity");
  dependencies_.insert(KeyTag{ por_key_, tag });

  depth_key_ = Keys::readKey(plist, domain, "depth", "depth");
  dependencies_.insert(KeyTag{ depth_key_, tag });

  cv_key_ = Keys::readKey(plist, domain, "cell volume", "cell_volume");
  dependencies_.insert(KeyTag{ cv_key_, tag });

  q10_ = plist_.get<double>("Q10 [-]", 2.0);
}


Teuchos::RCP<Evaluator>
CarbonDecomposeRateEvaluator::Clone() const
{
  return Teuchos::rcp(new CarbonDecomposeRateEvaluator(*this));
}


void
CarbonDecomposeRateEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  Epetra_MultiVector& res_c = *result[0]->ViewComponent("cell", false);

  const auto& temp_c = *S.Get<CompositeVector>(temp_key_, tag).ViewComponent("cell", false);
  const auto& pres_c = *S.Get<CompositeVector>(pres_key_, tag).ViewComponent("cell", false);
  const auto& por_c = *S.Get<CompositeVector>(por_key_, tag).ViewComponent("cell", false);
  const auto& depth_c = *S.Get<CompositeVector>(depth_key_, tag).ViewComponent("cell", false);
  const auto& cv_c = *S.Get<CompositeVector>(cv_key_, tag).ViewComponent("cell", false);

  for (int c = 0; c != res_c.MyLength(); ++c) {
    if (temp_c[0][c] >= 273.15) {
      double f_temp = Func_Temp(temp_c[0][c], q10_);
      double f_depth = Func_Depth(depth_c[0][c]);
      double f_pres_temp = Func_TempPres(temp_c[0][c], pres_c[0][c]);
      res_c[0][c] = (f_temp * f_depth * f_pres_temp) * (1 - por_c[0][c]) * cv_c[0][c];
    } else {
      res_c[0][c] = 0.;
    }
  }
}


double
CarbonDecomposeRateEvaluator::Func_Temp(double temp, double q10) const
{
  return std::pow(q10, (temp - 293.15) / 10.0);
}


double
CarbonDecomposeRateEvaluator::Func_Depth(double depth) const
{
  return depth < 1.0 ? 1 : std::exp(-0.5 * (depth - 1));
}


double
CarbonDecomposeRateEvaluator::Func_TempPres(double temp, double pres) const
{
  double p_min = -1.0e7;
  double p_max = -1.0e4;
  double p_atm = 101325.;

  double pn_star = pres - p_atm;
  double pn = std::min(0.0, std::max(pn_star, p_min));

  double f = -1.e18;
  if (pn >= p_max) {
    f = pn / p_max;
  } else {
    AMANZI_ASSERT(pn != 0.);
    f = std::log(p_min / pn) / std::log(p_min / p_max);
  }
  return f;
}

} // namespace Relations
} // namespace Amanzi
