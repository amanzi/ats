/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The elevation evaluator gets the subsurface temperature and computes moisture content
  over time.

  Authors: Ahmad Jan (jana@ornl.gov)
*/

#include "carbon_decomposition_rate_evaluator.hh"

namespace Amanzi {
namespace Flow {


CarbonDecomposeRateEvaluator::CarbonDecomposeRateEvaluator(Teuchos::ParameterList& plist)
    : EvaluatorSecondaryMonotypeCV(plist)
{
  Tag tag = my_keys_.front().second;
  domain_ = Keys::getDomain(my_keys_.front().first); //surface_column domain
  Key domain_ss = Keys::readDomainHint(plist_, domain_, "surface", "subsurface");

  temp_key_ =  Keys::readKey(plist, domain_ss, "temperature","temperature");
  dependencies_.insert(KeyTag{temp_key_, tag});

  pres_key_ =  Keys::readKey(plist, domain_ss, "pressure", "pressure");
  dependencies_.insert(KeyTag{pres_key_, tag});

  por_key_ =  Keys::readKey(plist, domain_ss, "porosity", "porosity");
  dependencies_.insert(KeyTag{por_key_, tag});

  cv_key_ =  Keys::readKey(plist, domain_ss, "cell volume","cell_volume");
  dependencies_.insert(KeyTag{cv_key_, tag});

  //trans_width_ =  plist_.get<double>("transition width [K]", 0.2);
  q10_ =  plist_.get<double>("Q10 [-]", 2.0);
}


Teuchos::RCP<Evaluator>
CarbonDecomposeRateEvaluator::Clone() const
{
  return Teuchos::rcp(new CarbonDecomposeRateEvaluator(*this));
}


void
CarbonDecomposeRateEvaluator::Evaluate_(const State& S,
        const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  Epetra_MultiVector& res_c = *result[0]->ViewComponent("cell",false);
  AMANZI_ASSERT(res_c.MyLength() == 0); // this PK only valid on column mesh

  const auto& temp_c = *S.Get<CompositeVector>(temp_key_, tag).ViewComponent("cell", false);
  const auto& pres_c = *S.Get<CompositeVector>(pres_key_, tag).ViewComponent("cell", false);
  const auto& por_c = *S.Get<CompositeVector>(por_key_, tag).ViewComponent("cell", false);
  const auto& vol_c = *S.Get<CompositeVector>(cv_key_, tag).ViewComponent("cell", false);

  std::string domain_ss = Keys::getDomain(temp_key_);
  const auto& top_z_centroid = S.GetMesh(domain_ss)->face_centroid(0);
  AmanziGeometry::Point z_up_centroid(top_z_centroid);
  AmanziGeometry::Point z_down_centroid(top_z_centroid);
  AmanziGeometry::Point z_depth(top_z_centroid);

  int col_cells = temp_c.MyLength();
  double col_sum = 0;

  for (int i=0; i!=col_cells; ++i) {
    if (temp_c[0][i] >= 273.15) {

      z_up_centroid = S.GetMesh(domain_ss)->face_centroid(i);
      z_down_centroid = S.GetMesh(domain_ss)->face_centroid(i+1);

      double dz = z_up_centroid[2] - z_down_centroid[2];
      double depth = top_z_centroid[2] - z_down_centroid[2];

      double f_temp = std::pow(q10_,(temp_c[0][i]-293.15)/10.0);
      double f_depth = depth < 1.0 ? 1 : std::exp(-0.5*(depth-1));
      double f_pres_temp = Func_TempPres(temp_c[0][i],pres_c[0][i]);

      double soil = (1-por_c[0][i]) * vol_c[0][i];
      col_sum += (f_temp * f_depth * f_pres_temp * dz) * soil;
    }
  }
  res_c[0][0] = col_sum;
}

void
CarbonDecomposeRateEvaluator::EvaluatePartialDerivative_(const State& S,
               const Key& wrt_key, const Tag& wrt_tag, const std::vector<CompositeVector*>& result)
{}


// Custom EnsureCompatibility forces this to be updated once.
bool
CarbonDecomposeRateEvaluator::Update(State& S, const Key& request)
{
  bool changed = EvaluatorSecondaryMonotypeCV::Update(S,request);

  if (!updated_once_) {
    Update_(S);
    updated_once_ = true;
    return true;
  }
  return changed;
}

void
CarbonDecomposeRateEvaluator::EnsureCompatibility(State& S)
{
  // note, no derivs are valid here
  EnsureCompatibility_ClaimOwnership_(S);
  EnsureCompatibility_Flags_(S);
  EnsureCompatibility_DepEvals_(S);

  CompositeVectorSpace fac;
  fac.AddComponent("cell", AmanziMesh::CELL, 1);
  EnsureCompatibility_DepsFromFac_(S, fac, true);

  EnsureCompatibility_DepEnsureCompatibility_(S);
}


double
CarbonDecomposeRateEvaluator::Func_TempPres(double temp, double pres)
{
  double p_min = -1.0e7;
  double p_max = -1.0e4;
  double p_atm = 101325.;

  double pn_star = pres - p_atm;
  double pn = std::min(0.0, std::max(pn_star,p_min));

  double f = -1.e18;
  if (pn >= p_max) {
    f = pn / p_max;
  } else {
    AMANZI_ASSERT(pn != 0.);
    f = std::log(p_min/pn) / std::log(p_min/p_max);
  }
  return f;
}

} //namespace
} //namespace
