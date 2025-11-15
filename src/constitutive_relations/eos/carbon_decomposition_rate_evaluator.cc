/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ahmad Jan (jana@ornl.gov)
           Bo Gao (gaob@ornl.gov)
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

  sat_liq_key_ = Keys::readKey(plist, domain_, "saturation liquid", "saturation_liquid");
  dependencies_.insert(KeyTag{ sat_liq_key_, tag });

  sat_gas_key_ = Keys::readKey(plist, domain_, "saturation gas", "saturation_gas");
  dependencies_.insert(KeyTag{ sat_gas_key_, tag });

  por_key_ = Keys::readKey(plist, domain_, "porosity", "porosity");
  dependencies_.insert(KeyTag{ por_key_, tag });

  depth_key_ = Keys::readKey(plist, domain_, "depth", "depth");
  dependencies_.insert(KeyTag{ depth_key_, tag });

  subsidence_key_ = Keys::readKey(plist, domain_surf_, "subsidence", "subsidence");
  dependencies_.insert(KeyTag{ subsidence_key_, tag });

  q10_ = plist_.get<double>("Q10 [-]", 2.0);

  is_func_temp_ = plist_.get<bool>("use temperature dependent coefficient", true);
  is_func_depth_ = plist_.get<bool>("use depth dependent coefficient", false);
  is_func_pres_ = plist_.get<bool>("use pressure dependent coefficient", false);
  is_func_liq_ = plist_.get<bool>("use liquid saturation dependent coefficient", true);
  is_func_gas_ = plist_.get<bool>("use gas saturation dependent coefficient", true);
  is_scaling_down_ = plist_.get<bool>("scale down with increase of moisture", false);
  is_threshold_temp_ = plist_.get<bool>("count above freezing", true);
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
  const auto& sat_liq_c = *S.Get<CompositeVector>(sat_liq_key_, tag).ViewComponent("cell", false);
  const auto& sat_gas_c = *S.Get<CompositeVector>(sat_gas_key_, tag).ViewComponent("cell", false);
  const auto& subsidence_c = *S.Get<CompositeVector>(subsidence_key_, tag).ViewComponent("cell", false);

  const auto& mesh = *S.GetMesh(domain_);
  const auto& mesh_surf = *S.GetMesh(domain_surf_);

  double threshold_temp = is_threshold_temp_ ? 273.15 : -1;
  for (int col = 0; col != mesh.columns.num_columns_owned; ++col) {
    const auto& col_cells = mesh.columns.getCells(col);
    for (int c = 0; c != col_cells.size(); ++c) {
      if (temp_c[0][col_cells[c]] >= threshold_temp) {
        double dz = mesh.getCellVolume(col_cells[c]) / mesh_surf.getCellVolume(col);
        double f_temp = is_func_temp_ ? Func_Temp(temp_c[0][col_cells[c]], q10_) : 1.;
        double f_depth = is_func_depth_ ? Func_Depth(depth_c[0][col_cells[c]] + subsidence_c[0][col]) : 1.;
        double f_pres = is_func_pres_ ? Func_Pres(pres_c[0][col_cells[c]]) : 1.;
        double f_liq = is_func_liq_ ? Func_Satliq(sat_liq_c[0][col_cells[c]]) : 1.;
        double f_gas = is_func_gas_ ? Func_Satgas(sat_gas_c[0][col_cells[c]]) : 1.;
        res_c[0][col_cells[c]] = (f_temp * f_depth * f_pres * f_liq * f_gas * dz) * (1 - por_c[0][col_cells[c]]);
      } else {
        res_c[0][col_cells[c]] = 0.;
      }
    }
  }
}


void
CarbonDecomposeRateEvaluator::EnsureCompatibility_ToDeps_(State& S)
{
  const auto& fac = S.Require<CompositeVector, CompositeVectorSpace>(my_keys_.front().first,
                                                                     my_keys_.front().second);
  if (fac.Mesh() != Teuchos::null) {
    for (const auto& dep : dependencies_) {
      if (Keys::getDomain(dep.first) == Keys::getDomain(my_keys_.front().first)) {
        S.Require<CompositeVector, CompositeVectorSpace>(dep.first, dep.second).Update(fac);
      } else {
        CompositeVectorSpace dep_fac;
        dep_fac.SetMesh(S.GetMesh(Keys::getDomain(dep.first)))
          ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
        S.Require<CompositeVector, CompositeVectorSpace>(dep.first, dep.second).Update(dep_fac);
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
CarbonDecomposeRateEvaluator::Func_Pres(double pres) const
{
  double p_min = -2.5e6;
  double p_max = -2.0e3;
  double p_atm = 101325.;

  double pn_star = pres - p_atm;
  double pn = std::min(0.0, std::max(pn_star, p_min));

  double f = -1.e18;
  if (pn >= p_max) {
    if (is_scaling_down_) {
      f = pn / p_max;
    } else {
      f = 1.;
    }
  } else {
    AMANZI_ASSERT(pn != 0.);
    f = std::log(p_min / pn) / std::log(p_min / p_max);
  }
  return f;
}


double
CarbonDecomposeRateEvaluator::Func_Satliq(double sat) const
{
  if (sat <= 0.7) {
    return std::pow(sat / 0.7, 1.5);
  } else {
    return std::pow((1 - sat) / 0.3, 1.5);
  }
}


double
CarbonDecomposeRateEvaluator::Func_Satgas(double sat) const
{
  if (sat <= 0.3) {
    return std::pow(sat / 0.3, 1.5);
  } else {
    return std::pow((1 - sat) / 0.7, 1.5);
  }
}

} // namespace Relations
} // namespace Amanzi
