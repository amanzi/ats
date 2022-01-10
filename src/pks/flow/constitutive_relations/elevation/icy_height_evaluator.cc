/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluator for determining height( rho, head )

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "icy_height_model.hh"
#include "icy_height_evaluator.hh"


namespace Amanzi {
namespace Flow {


IcyHeightEvaluator::IcyHeightEvaluator(Teuchos::ParameterList& plist) :
    HeightEvaluator(plist) {

  Key domain = Keys::getDomain(my_key_);
  // my extra dependencies
  dens_ice_key_ = plist_.get<std::string>("ice mass density key", Keys::getKey(domain,"mass_density_ice"));  
  dependencies_.insert(dens_ice_key_);

  unfrozen_frac_key_ = plist_.get<std::string>("unfrozen fraction key", Keys::getKey(domain,"unfrozen_fraction"));

  dependencies_.insert(unfrozen_frac_key_);

  // model
  Teuchos::ParameterList model_plist = plist_.sublist("height model parameters");
  icy_model_ = Teuchos::rcp(new IcyHeightModel(model_plist));

}


IcyHeightEvaluator::IcyHeightEvaluator(const IcyHeightEvaluator& other) :
    HeightEvaluator(other),
    dens_ice_key_(other.dens_ice_key_),
    unfrozen_frac_key_(other.unfrozen_frac_key_),
    icy_model_(other.icy_model_) {}


Teuchos::RCP<FieldEvaluator>
IcyHeightEvaluator::Clone() const {
  return Teuchos::rcp(new IcyHeightEvaluator(*this));
}


void IcyHeightEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {

  Teuchos::RCP<const CompositeVector> pres = S->GetPtr<CompositeVector>(pres_key_);

  // this is rather hacky.  surface_pressure is a mixed field vector -- it has
  // pressure on cells and ponded depth on faces.
  // -- copy the faces over directly
  if (result->HasComponent("face"))
    *result->ViewComponent("face",false) = *pres->ViewComponent("face",false);

  // -- cells need the function eval
  const Epetra_MultiVector& res_c = *result->ViewComponent("cell",false);
  const Epetra_MultiVector& pres_c = *pres->ViewComponent("cell",false);
  const Epetra_MultiVector& rho_l = *S->GetPtr<CompositeVector>(dens_key_)
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& rho_i = *S->GetPtr<CompositeVector>(dens_ice_key_)
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& eta = *S->GetPtr<CompositeVector>(unfrozen_frac_key_)
      ->ViewComponent("cell",false);

  const double& p_atm = *S->GetScalarData(patm_key_);
  const Epetra_Vector& gravity = *S->GetConstantVectorData(gravity_key_);
  double gz = -gravity[2];

  int ncells = res_c.MyLength();
  if (bar_) {
    for (int c=0; c!=ncells; ++c) {
      res_c[0][c] = icy_model_->Height(pres_c[0][c], eta[0][c],
              rho_l[0][c], rho_i[0][c], p_atm, gz);
    }
  } else {
    for (int c=0; c!=ncells; ++c) {
      res_c[0][c] = pres_c[0][c] < p_atm ? 0. :
          icy_model_->Height(pres_c[0][c], eta[0][c],
                             rho_l[0][c], rho_i[0][c], p_atm, gz);
    }
  }
}


void IcyHeightEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
  // this is rather hacky.  surface_pressure is a mixed field vector -- it has
  // pressure on cells and ponded depth on faces.
  // -- NO FACE DERIVATIVES
  //  result->ViewComponent("face",false)->PutScalar(1.0);

  // -- cells need the function eval
  const Epetra_MultiVector& res_c = *result->ViewComponent("cell",false);
  const Epetra_MultiVector& pres_c = *S->GetPtr<CompositeVector>(pres_key_)
      ->ViewComponent("cell",false);
   const Epetra_MultiVector& rho_l = *S->GetPtr<CompositeVector>(dens_key_)
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& rho_i = *S->GetPtr<CompositeVector>(dens_ice_key_)
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& eta = *S->GetPtr<CompositeVector>(unfrozen_frac_key_)
      ->ViewComponent("cell",false);

  const double& p_atm = *S->GetScalarData(patm_key_);
  const Epetra_Vector& gravity = *S->GetConstantVectorData(gravity_key_);
  double gz = -gravity[2];  // check this

  // For derivatives, the height is always assumed to be non-negative.  If it
  // is negative, the term gets zeroed later.
  if (bar_) {
    if (wrt_key == pres_key_) {
      int ncells = res_c.MyLength();
      for (int c=0; c!=ncells; ++c) {
        res_c[0][c] = icy_model_->DHeightDPressure(pres_c[0][c], eta[0][c],
                rho_l[0][c], rho_i[0][c], p_atm, gz);
      }
    } else if (wrt_key == dens_key_) {
      int ncells = res_c.MyLength();
      for (int c=0; c!=ncells; ++c) {
        res_c[0][c] = icy_model_->DHeightDRho_l(pres_c[0][c], eta[0][c],
                rho_l[0][c], rho_i[0][c], p_atm, gz);
      }
    } else if (wrt_key == dens_ice_key_) {
      int ncells = res_c.MyLength();
      for (int c=0; c!=ncells; ++c) {
        res_c[0][c] = icy_model_->DHeightDRho_i(pres_c[0][c], eta[0][c],
                rho_l[0][c], rho_i[0][c], p_atm, gz);
      }
    } else if (wrt_key == unfrozen_frac_key_) {
      int ncells = res_c.MyLength();
      for (int c=0; c!=ncells; ++c) {
        res_c[0][c] = icy_model_->DHeightDEta(pres_c[0][c], eta[0][c],
                rho_l[0][c], rho_i[0][c], p_atm, gz);
      }
    } else {
      AMANZI_ASSERT(0);
    }
  } else {
    if (wrt_key == pres_key_) {
      int ncells = res_c.MyLength();
      for (int c=0; c!=ncells; ++c) {
        res_c[0][c] =  pres_c[0][c] < p_atm ? 0. :
            icy_model_->DHeightDPressure(pres_c[0][c], eta[0][c],
                rho_l[0][c], rho_i[0][c], p_atm, gz);
      }
    } else if (wrt_key == dens_key_) {
      int ncells = res_c.MyLength();
      for (int c=0; c!=ncells; ++c) {
        res_c[0][c] =  pres_c[0][c] < p_atm ? 0. :
            icy_model_->DHeightDRho_l(pres_c[0][c], eta[0][c],
                rho_l[0][c], rho_i[0][c], p_atm, gz);
      }
    } else if (wrt_key == dens_ice_key_) {
      int ncells = res_c.MyLength();
      for (int c=0; c!=ncells; ++c) {
        res_c[0][c] =  pres_c[0][c] < p_atm ? 0. :
            icy_model_->DHeightDRho_i(pres_c[0][c], eta[0][c],
                rho_l[0][c], rho_i[0][c], p_atm, gz);
      }
    } else if (wrt_key == unfrozen_frac_key_) {
      int ncells = res_c.MyLength();
      for (int c=0; c!=ncells; ++c) {
        res_c[0][c] =  pres_c[0][c] < p_atm ? 0. :
            icy_model_->DHeightDEta(pres_c[0][c], eta[0][c],
                rho_l[0][c], rho_i[0][c], p_atm, gz);
      }
    } else {
      AMANZI_ASSERT(0);
    }
  }
}



} //namespace
} //namespace
