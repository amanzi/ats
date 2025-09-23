/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatsky(dasvyat@lanl.gov)
*/

/*
  Evaluator for determining amount of water on the surface based on density and pressure

*/

#include "overland_pressure_multicomponent_water_content_evaluator.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {


OverlandPressureMulticomponentWaterContentEvaluator::
  OverlandPressureMulticomponentWaterContentEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  bar_ = plist_.get<bool>("allow negative water content", false);
  rollover_ = plist_.get<double>("water content rollover", 0.);

  Key domain_name = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  // my dependencies
  pres_key_ = Keys::readKey(plist_, domain_name, "pressure", "pressure");
  dependencies_.insert(KeyTag{ pres_key_, tag });

  mass_dens_key_ = Keys::readKey(plist_, domain_name, "mass density", "mass_density_liquid");
  dependencies_.insert(KeyTag{ mass_dens_key_, tag });
  molar_dens_key_ = Keys::readKey(plist_, domain_name, "molar density", "molar_density_liquid");
  dependencies_.insert(KeyTag{ molar_dens_key_, tag });

  cv_key_ = Keys::readKey(plist_, domain_name, "cell volume", "cell_volume");
  dependencies_.insert(KeyTag{ cv_key_, tag });
}


Teuchos::RCP<Evaluator>
OverlandPressureMulticomponentWaterContentEvaluator::Clone() const
{
  return Teuchos::rcp(new OverlandPressureMulticomponentWaterContentEvaluator(*this));
}


void
OverlandPressureMulticomponentWaterContentEvaluator::Evaluate_(
  const State& S,
  const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  Epetra_MultiVector& res = *result[0]->ViewComponent("cell", false);
  const Epetra_MultiVector& pres =
    *S.GetPtr<CompositeVector>(pres_key_, tag)->ViewComponent("cell", false);
  const Epetra_MultiVector& molar_dens =
    *S.GetPtr<CompositeVector>(molar_dens_key_, tag)->ViewComponent("cell", false);
  const Epetra_MultiVector& mass_dens =
    *S.GetPtr<CompositeVector>(mass_dens_key_, tag)->ViewComponent("cell", false);


  const Epetra_MultiVector& cv =
    *S.GetPtr<CompositeVector>(cv_key_, tag)->ViewComponent("cell", false);

  const double& p_atm = S.Get<double>("atmospheric_pressure", Tags::DEFAULT);
  const auto& gravity = S.Get<AmanziGeometry::Point>("gravity", Tags::DEFAULT);
  double gz = -gravity[2]; // check this

  int ncells = res.MyLength();
  if (bar_) {
    for (int c = 0; c != ncells; ++c) {
      res[0][c] = cv[0][c] * molar_dens[0][c] * (pres[0][c] - p_atm) / (mass_dens[0][c] * gz);
    }
  } else if (rollover_ > 0.) {
    for (int c = 0; c != ncells; ++c) {
      double dp = pres[0][c] - p_atm;
      double dp_eff = dp < 0.        ? 0. :
                      dp < rollover_ ? dp * dp / (2 * rollover_) :
                                       dp - rollover_ / 2.;
      res[0][c] = cv[0][c] * molar_dens[0][c] * dp_eff / (mass_dens[0][c] * gz);
    }
  } else {
    for (int c = 0; c != ncells; ++c) {
      res[0][c] = pres[0][c] < p_atm ?
                    0. :
                    cv[0][c] * molar_dens[0][c] * (pres[0][c] - p_atm) / (mass_dens[0][c] * gz);
    }
  }
}


void
OverlandPressureMulticomponentWaterContentEvaluator::EvaluatePartialDerivative_(
  const State& S,
  const Key& wrt_key,
  const Tag& wrt_tag,
  const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  AMANZI_ASSERT(wrt_key == pres_key_);

  Epetra_MultiVector& res = *result[0]->ViewComponent("cell", false);
  const Epetra_MultiVector& pres =
    *S.GetPtr<CompositeVector>(pres_key_, tag)->ViewComponent("cell", false);
  const Epetra_MultiVector& molar_dens =
    *S.GetPtr<CompositeVector>(molar_dens_key_, tag)->ViewComponent("cell", false);
  const Epetra_MultiVector& mass_dens =
    *S.GetPtr<CompositeVector>(mass_dens_key_, tag)->ViewComponent("cell", false);

  const Epetra_MultiVector& cv =
    *S.GetPtr<CompositeVector>(cv_key_, tag)->ViewComponent("cell", false);

  const double& p_atm = S.Get<double>("atmospheric_pressure", Tags::DEFAULT);
  const auto& gravity = S.Get<AmanziGeometry::Point>("gravity", Tags::DEFAULT);
  double gz = -gravity[2]; // check this

  if (wrt_key == pres_key_) {
    int ncells = res.MyLength();
    if (bar_) {
      for (int c = 0; c != ncells; ++c) {
        res[0][c] = cv[0][c] * molar_dens[0][c] / (mass_dens[0][c] * gz);
      }
    } else if (rollover_ > 0.) {
      for (int c = 0; c != ncells; ++c) {
        double dp = pres[0][c] - p_atm;
        double ddp_eff = dp < 0. ? 0. : dp < rollover_ ? dp / rollover_ : 1.;
        res[0][c] = cv[0][c] * ddp_eff * molar_dens[0][c] / (mass_dens[0][c] * gz);
      }
    } else {
      for (int c = 0; c != ncells; ++c) {
        res[0][c] = pres[0][c] < p_atm ? 0. : cv[0][c] * molar_dens[0][c] / (mass_dens[0][c] * gz);
      }
    }
  } else {
    res.PutScalar(0.);
  }
}


} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi
