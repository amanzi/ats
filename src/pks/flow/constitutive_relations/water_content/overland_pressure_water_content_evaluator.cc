/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluator for determining height( rho, head )

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "overland_pressure_water_content_evaluator.hh"

namespace Amanzi {
namespace Flow {


OverlandPressureWaterContentEvaluator::OverlandPressureWaterContentEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {
  M_ = plist_.get<double>("molar mass", 0.0180153);
  bar_ = plist_.get<bool>("allow negative water content", false);
  rollover_ = plist_.get<double>("water content rollover", 0.);

  Key domain = Keys::getDomain(my_key_);

  // my dependencies
  pres_key_ = plist_.get<std::string>("pressure key", Keys::getKey(domain, "pressure"));
  dependencies_.insert(pres_key_);
  cv_key_ = plist_.get<std::string>("cell volume key", Keys::getKey(domain, "cell_volume"));
  dependencies_.insert(cv_key_);
}


OverlandPressureWaterContentEvaluator::OverlandPressureWaterContentEvaluator(const OverlandPressureWaterContentEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    pres_key_(other.pres_key_),
    M_(other.M_),
    bar_(other.bar_),
    rollover_(other.rollover_),
    cv_key_(other.cv_key_)
{}


Teuchos::RCP<FieldEvaluator>
OverlandPressureWaterContentEvaluator::Clone() const {
  return Teuchos::rcp(new OverlandPressureWaterContentEvaluator(*this));
}


void OverlandPressureWaterContentEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {

  Epetra_MultiVector& res = *result->ViewComponent("cell",false);
  const Epetra_MultiVector& pres = *S->GetPtr<CompositeVector>(pres_key_)
      ->ViewComponent("cell",false);

  const Epetra_MultiVector& cv = *S->GetPtr<CompositeVector>(cv_key_)
      ->ViewComponent("cell",false);

  const double& p_atm = *S->GetScalarData("atmospheric_pressure");
  const Epetra_Vector& gravity = *S->GetConstantVectorData("gravity");
  double gz = -gravity[2];  // check this

  int ncells = res.MyLength();
  if (bar_) {
    for (int c=0; c!=ncells; ++c) {
      res[0][c] = cv[0][c] * (pres[0][c] - p_atm) / (gz * M_);
    }
  } else if (rollover_ > 0.) {
    for (int c=0; c!=ncells; ++c) {
      double dp = pres[0][c] - p_atm;
      double dp_eff = dp < 0. ? 0. :
          dp < rollover_ ?
            dp*dp/(2*rollover_) :
            dp - rollover_/2.;
      res[0][c] = cv[0][c] * dp_eff / (gz * M_);
    }
  } else {
    for (int c=0; c!=ncells; ++c) {
      res[0][c] = pres[0][c] < p_atm ? 0. :
          cv[0][c] * (pres[0][c] - p_atm) / (gz * M_);
    }
  }
}


void OverlandPressureWaterContentEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
  AMANZI_ASSERT(wrt_key == pres_key_);

  Epetra_MultiVector& res = *result->ViewComponent("cell",false);
  const Epetra_MultiVector& pres = *S->GetPtr<CompositeVector>(pres_key_)
      ->ViewComponent("cell",false);

  const Epetra_MultiVector& cv = *S->GetPtr<CompositeVector>(cv_key_)
      ->ViewComponent("cell",false);

  const double& p_atm = *S->GetScalarData("atmospheric_pressure");
  const Epetra_Vector& gravity = *S->GetConstantVectorData("gravity");
  double gz = -gravity[2];  // check this

  int ncells = res.MyLength();
  if (bar_) {
    for (int c=0; c!=ncells; ++c) {
      res[0][c] = cv[0][c] / (gz * M_);
    }
  } else if (rollover_ > 0.) {
    for (int c=0; c!=ncells; ++c) {
      double dp = pres[0][c] - p_atm;
      double ddp_eff = dp < 0. ? 0. :
          dp < rollover_ ? dp/rollover_ : 1.;
      res[0][c] = cv[0][c] * ddp_eff / (gz * M_);
    }
  } else {
    for (int c=0; c!=ncells; ++c) {
      res[0][c] = pres[0][c] < p_atm ? 0. :
          cv[0][c] / (gz * M_);
    }
  }
}


} //namespace
} //namespace
