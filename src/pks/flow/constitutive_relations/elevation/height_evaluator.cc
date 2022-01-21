/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluator for determining height( rho, head )

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "height_model.hh"
#include "height_evaluator.hh"


namespace Amanzi {
namespace Flow {


HeightEvaluator::HeightEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  bar_ = plist_.get<bool>("allow negative ponded depth", false);
  Key domain = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  // my dependencies
  dens_key_ = Keys::readKey(plist_, domain, "mass density", "mass_density_liquid");
  dependencies_.insert(KeyTag{dens_key_, tag});

  pres_key_ = Keys::readKey(plist_, domain, "pressure", "pressure");
  dependencies_.insert(KeyTag{pres_key_, tag});

  // model
  Teuchos::ParameterList model_plist = plist_.sublist("height model parameters");
  model_ = Teuchos::rcp(new HeightModel(model_plist));
}


Teuchos::RCP<Evaluator>
HeightEvaluator::Clone() const {
  return Teuchos::rcp(new HeightEvaluator(*this));
}


void HeightEvaluator::EnsureCompatibility(State& S)
{
  EnsureCompatibility_ClaimOwnership_(S);
  EnsureCompatibility_Flags_(S);
  EnsureCompatibility_Derivs_(S);
  EnsureCompatibility_DepEvals_(S);

  auto akeytag = my_keys_[0];
  const auto& my_fac = S.Require<CompositeVector,CompositeVectorSpace>(akeytag.first, akeytag.second);
  if (my_fac.Mesh() != Teuchos::null) {
    // Create an unowned factory to check my dependencies.  This is done
    // manually here because we do NOT want faces, despite having faces in
    // my_key.  The faces will get updated directly from the mixed field.
    CompositeVectorSpace dep_fac;
    dep_fac.SetOwned(false);
    dep_fac.SetGhosted(my_fac.Ghosted());
    dep_fac.SetMesh(my_fac.Mesh());
    dep_fac.AddComponent("cell", AmanziMesh::CELL, 1);

    EnsureCompatibility_DepsFromFac_(S, dep_fac);
  }

  EnsureCompatibility_DepDerivs_(S);
  EnsureCompatibility_DepEnsureCompatibility_(S);
}


void HeightEvaluator::Evaluate_(const State& S,
        const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> pres = S.GetPtr<CompositeVector>(pres_key_, tag);

  // this is rather hacky.  surface_pressure is a mixed field vector -- it has
  // pressure on cells and ponded depth on faces.
  // -- copy the faces over directly
  if (result[0]->HasComponent("face") && pres->HasComponent("face"))
    *result[0]->ViewComponent("face",false) = *pres->ViewComponent("face",false);

  // -- cells need the function eval
  const Epetra_MultiVector& res_c = *result[0]->ViewComponent("cell",false);
  const Epetra_MultiVector& pres_c = *pres->ViewComponent("cell",false);
  const Epetra_MultiVector& rho = *S.GetPtr<CompositeVector>(dens_key_, tag)
      ->ViewComponent("cell",false);

  double p_atm = S.Get<double>("atmospheric_pressure");
  const AmanziGeometry::Point& gravity = S.Get<AmanziGeometry::Point>("gravity");
  double gz = -gravity[2];

  int ncells = res_c.MyLength();
  if (bar_) {
    for (int c=0; c!=ncells; ++c) {
      res_c[0][c] = model_->Height(pres_c[0][c], rho[0][c], p_atm, gz);
    }
  } else {
    for (int c=0; c!=ncells; ++c) {
      res_c[0][c] = pres_c[0][c] < p_atm ? 0. :
          model_->Height(pres_c[0][c], rho[0][c], p_atm, gz);
    }
  }
}


void HeightEvaluator::EvaluatePartialDerivative_(const State& S,
        const Key& wrt_key, const Tag& wrt_tag,
        const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;

  // -- cells need the function eval
  const Epetra_MultiVector& res_c = *result[0]->ViewComponent("cell",false);
  const Epetra_MultiVector& pres_c = *S.GetPtr<CompositeVector>(pres_key_, tag)
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& rho = *S.GetPtr<CompositeVector>(dens_key_, tag)
      ->ViewComponent("cell",false);

  double p_atm = S.Get<double>("atmospheric_pressure");
  const AmanziGeometry::Point& gravity = S.Get<AmanziGeometry::Point>("gravity");
  double gz = -gravity[2];

  if (wrt_key == pres_key_) {
    int ncells = res_c.MyLength();
    if (bar_) {
      for (int c=0; c!=ncells; ++c) {
        res_c[0][c] = model_->DHeightDPressure(pres_c[0][c], rho[0][c], p_atm, gz);
      }
    } else {
      for (int c=0; c!=ncells; ++c) {
        res_c[0][c] =  pres_c[0][c] < p_atm ? 0. :
            model_->DHeightDPressure(pres_c[0][c], rho[0][c], p_atm, gz);
      }
    }
  } else if (wrt_key == dens_key_) {
    int ncells = res_c.MyLength();
    if (bar_) {
      for (int c=0; c!=ncells; ++c) {
        res_c[0][c] = model_->DHeightDRho(pres_c[0][c], rho[0][c], p_atm, gz);
      }
    } else {
      for (int c=0; c!=ncells; ++c) {
        res_c[0][c] =  pres_c[0][c] < p_atm ? 0. :
            model_->DHeightDRho(pres_c[0][c], rho[0][c], p_atm, gz);
      }
    }
  } else {
    AMANZI_ASSERT(0);
  }
}



} //namespace
} //namespace
