/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "icy_height_model.hh"
#include "icy_height_evaluator.hh"


namespace Amanzi {
namespace ATS_Physics {
namespace Flow {


IcyHeightEvaluator::IcyHeightEvaluator(Teuchos::ParameterList& plist)
  : HeightEvaluator(plist)
{
  Key domain = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  // my extra dependencies
  dens_ice_key_ = Keys::readKey(plist_, domain, "ice mass density", "mass_density_ice");
  dependencies_.insert(KeyTag{ dens_ice_key_, tag });

  unfrozen_frac_key_ = Keys::readKey(plist_, domain, "unfrozen fraction", "unfrozen_fraction");
  dependencies_.insert(KeyTag{ unfrozen_frac_key_, tag });

  // model
  Teuchos::ParameterList model_plist = plist_.sublist("height model parameters");
  icy_model_ = Teuchos::rcp(new IcyHeightModel(model_plist));
}


Teuchos::RCP<Evaluator>
IcyHeightEvaluator::Clone() const
{
  return Teuchos::rcp(new IcyHeightEvaluator(*this));
}


void
IcyHeightEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> pres = S.GetPtr<CompositeVector>(pres_key_, tag);

  // this is rather hacky.  surface_pressure is a mixed field vector -- it has
  // pressure on cells and ponded depth on faces.
  // -- copy the faces over directly
  if (result[0]->HasComponent("face"))
    *result[0]->ViewComponent("face", false) = *pres->ViewComponent("face", false);

  // -- cells need the function eval
  const Epetra_MultiVector& res_c = *result[0]->ViewComponent("cell", false);
  const Epetra_MultiVector& pres_c = *pres->ViewComponent("cell", false);
  const Epetra_MultiVector& rho_l =
    *S.GetPtr<CompositeVector>(dens_key_, tag)->ViewComponent("cell", false);
  const Epetra_MultiVector& rho_i =
    *S.GetPtr<CompositeVector>(dens_ice_key_, tag)->ViewComponent("cell", false);
  const Epetra_MultiVector& eta =
    *S.GetPtr<CompositeVector>(unfrozen_frac_key_, tag)->ViewComponent("cell", false);

  double p_atm = S.Get<double>("atmospheric_pressure", Tags::DEFAULT);
  const AmanziGeometry::Point& gravity = S.Get<AmanziGeometry::Point>("gravity", Tags::DEFAULT);
  double gz = -gravity[2];

  int ncells = res_c.MyLength();
  if (bar_) {
    for (int c = 0; c != ncells; ++c) {
      res_c[0][c] =
        icy_model_->Height(pres_c[0][c], eta[0][c], rho_l[0][c], rho_i[0][c], p_atm, gz);
    }
  } else {
    for (int c = 0; c != ncells; ++c) {
      res_c[0][c] =
        pres_c[0][c] < p_atm ?
          0. :
          icy_model_->Height(pres_c[0][c], eta[0][c], rho_l[0][c], rho_i[0][c], p_atm, gz);
    }
  }
}


void
IcyHeightEvaluator::EvaluatePartialDerivative_(const State& S,
                                               const Key& wrt_key,
                                               const Tag& wrt_tag,
                                               const std::vector<CompositeVector*>& result)
{
  // this is rather hacky.  surface_pressure is a mixed field vector -- it has
  // pressure on cells and ponded depth on faces.
  // -- NO FACE DERIVATIVES
  Tag tag = my_keys_.front().second;

  // -- cells need the function eval
  const Epetra_MultiVector& res_c = *result[0]->ViewComponent("cell", false);
  const Epetra_MultiVector& pres_c =
    *S.GetPtr<CompositeVector>(pres_key_, tag)->ViewComponent("cell", false);
  const Epetra_MultiVector& rho_l =
    *S.GetPtr<CompositeVector>(dens_key_, tag)->ViewComponent("cell", false);
  const Epetra_MultiVector& rho_i =
    *S.GetPtr<CompositeVector>(dens_ice_key_, tag)->ViewComponent("cell", false);
  const Epetra_MultiVector& eta =
    *S.GetPtr<CompositeVector>(unfrozen_frac_key_, tag)->ViewComponent("cell", false);

  double p_atm = S.Get<double>("atmospheric_pressure", Tags::DEFAULT);
  const AmanziGeometry::Point& gravity = S.Get<AmanziGeometry::Point>("gravity", Tags::DEFAULT);
  double gz = -gravity[2];

  // For derivatives, the height is always assumed to be non-negative.  If it
  // is negative, the term gets zeroed later.
  if (bar_) {
    if (wrt_key == pres_key_) {
      int ncells = res_c.MyLength();
      for (int c = 0; c != ncells; ++c) {
        res_c[0][c] = icy_model_->DHeightDPressure(
          pres_c[0][c], eta[0][c], rho_l[0][c], rho_i[0][c], p_atm, gz);
      }
    } else if (wrt_key == dens_key_) {
      int ncells = res_c.MyLength();
      for (int c = 0; c != ncells; ++c) {
        res_c[0][c] =
          icy_model_->DHeightDRho_l(pres_c[0][c], eta[0][c], rho_l[0][c], rho_i[0][c], p_atm, gz);
      }
    } else if (wrt_key == dens_ice_key_) {
      int ncells = res_c.MyLength();
      for (int c = 0; c != ncells; ++c) {
        res_c[0][c] =
          icy_model_->DHeightDRho_i(pres_c[0][c], eta[0][c], rho_l[0][c], rho_i[0][c], p_atm, gz);
      }
    } else if (wrt_key == unfrozen_frac_key_) {
      int ncells = res_c.MyLength();
      for (int c = 0; c != ncells; ++c) {
        res_c[0][c] =
          icy_model_->DHeightDEta(pres_c[0][c], eta[0][c], rho_l[0][c], rho_i[0][c], p_atm, gz);
      }
    } else {
      AMANZI_ASSERT(0);
    }
  } else {
    if (wrt_key == pres_key_) {
      int ncells = res_c.MyLength();
      for (int c = 0; c != ncells; ++c) {
        res_c[0][c] = pres_c[0][c] < p_atm ?
                        0. :
                        icy_model_->DHeightDPressure(
                          pres_c[0][c], eta[0][c], rho_l[0][c], rho_i[0][c], p_atm, gz);
      }
    } else if (wrt_key == dens_key_) {
      int ncells = res_c.MyLength();
      for (int c = 0; c != ncells; ++c) {
        res_c[0][c] =
          pres_c[0][c] < p_atm ?
            0. :
            icy_model_->DHeightDRho_l(pres_c[0][c], eta[0][c], rho_l[0][c], rho_i[0][c], p_atm, gz);
      }
    } else if (wrt_key == dens_ice_key_) {
      int ncells = res_c.MyLength();
      for (int c = 0; c != ncells; ++c) {
        res_c[0][c] =
          pres_c[0][c] < p_atm ?
            0. :
            icy_model_->DHeightDRho_i(pres_c[0][c], eta[0][c], rho_l[0][c], rho_i[0][c], p_atm, gz);
      }
    } else if (wrt_key == unfrozen_frac_key_) {
      int ncells = res_c.MyLength();
      for (int c = 0; c != ncells; ++c) {
        res_c[0][c] =
          pres_c[0][c] < p_atm ?
            0. :
            icy_model_->DHeightDEta(pres_c[0][c], eta[0][c], rho_l[0][c], rho_i[0][c], p_atm, gz);
      }
    } else {
      AMANZI_ASSERT(0);
    }
  }
}


} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi
