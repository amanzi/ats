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
  domain_ = Keys::getDomain(my_keys_.front().first); // column, domain
  domain_surf_ = Keys::readDomainHint(plist, domain_, "subsurface", "surface");

  temp_key_ = Keys::readKey(plist, domain_, "temperature", "temperature");
  dependencies_.insert(KeyTag{ temp_key_, tag });

  pres_key_ = Keys::readKey(plist, domain_, "pressure", "pressure");
  dependencies_.insert(KeyTag{ pres_key_, tag });

  por_key_ = Keys::readKey(plist, domain_, "porosity", "porosity");
  dependencies_.insert(KeyTag{ por_key_, tag });

  depth_key_ = Keys::readKey(plist, domain_, "depth", "depth");
  dependencies_.insert(KeyTag{ depth_key_, tag });

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

  const auto& mesh = *S.GetMesh(domain_);
  const auto& mesh_surf = *S.GetMesh(domain_surf_);

  for (int col = 0; col != mesh.columns.num_columns_owned; ++col) {
    const auto& col_cells = mesh.columns.getCells(col);
    for (int c = 0; c != col_cells.size(); ++c) {
      if (temp_c[0][col_cells[c]] >= 273.15) {
        double dz = mesh.getCellVolume(col_cells[c]) / mesh_surf.getCellVolume(col);
        double f_temp = Func_Temp(temp_c[0][col_cells[c]], q10_);
        double f_depth = Func_Depth(depth_c[0][col_cells[c]]);
        double f_pres_temp = Func_TempPres(temp_c[0][col_cells[c]], pres_c[0][col_cells[c]]);
        res_c[0][col_cells[c]] =
          (f_temp * f_depth * f_pres_temp * dz) * (1 - por_c[0][col_cells[c]]);
      } else {
        res_c[0][col_cells[c]] = 0.;
      }
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
