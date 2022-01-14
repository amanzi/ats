/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluator for determining height( rho, head )

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "height_model.hh"
#include "height_evaluator.hh"


namespace Amanzi {
namespace Flow {


HeightEvaluator::HeightEvaluator(Teuchos::ParameterList& plist) :
    EvaluatorSecondaryMonotypeCV(plist) {
  bar_ = plist_.get<bool>("allow negative ponded depth", false);
  Key domain = Keys::getDomain(my_key_);

  // my dependencies
  dens_key_ = plist_.get<std::string>("mass density key", Keys::getKey(domain, "mass_density_liquid"));
  dependencies_.insert(dens_key_);

  pres_key_ = plist_.get<std::string>("pressure key", Keys::getKey(domain, "pressure"));
  dependencies_.insert(pres_key_);

  gravity_key_ = plist_.get<std::string>("gravity key", "gravity");
  patm_key_ = plist_.get<std::string>("atmospheric pressure key", "atmospheric_pressure");

  // model
  Teuchos::ParameterList model_plist = plist_.sublist("height model parameters");
  model_ = Teuchos::rcp(new HeightModel(model_plist));
}


HeightEvaluator::HeightEvaluator(const HeightEvaluator& other) :
    EvaluatorSecondaryMonotypeCV(other),
    dens_key_(other.dens_key_),
    pres_key_(other.pres_key_),
    gravity_key_(other.gravity_key_),
    patm_key_(other.patm_key_),
    model_(other.model_),
    bar_(other.bar_) {}


Teuchos::RCP<Evaluator>
HeightEvaluator::Clone() const {
  return Teuchos::rcp(new HeightEvaluator(*this));
}


void HeightEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S) {
  // Ensure my field exists.  Requirements should be already set.
  AMANZI_ASSERT(my_key_ != std::string(""));
  S->RequireGravity();
  Teuchos::RCP<CompositeVectorSpace> my_fac = S->Require<CompositeVector,CompositeVectorSpace>(my_key_, Tags::NEXT,  my_key_);
  my_fac->SetOwned(false);

  // check plist for vis or checkpointing control
  bool io_my_key = plist_.get<bool>(std::string("visualize ")+my_key_, true);
  S->GetField(my_key_, my_key_)->set_io_vis(io_my_key);
  bool checkpoint_my_key = plist_.get<bool>(std::string("checkpoint ")+my_key_, true);
  S->GetField(my_key_, my_key_)->set_io_checkpoint(checkpoint_my_key);

  // If my requirements have not yet been set, we'll have to hope they
  // get set by someone later.  For now just defer.
  if (my_fac->Mesh() != Teuchos::null) {
    // Create an unowned factory to check my dependencies.  This is done
    // manually here because we do NOT want faces, despite having faces in
    // my_key.  The faces will get updated directly from the mixed field.
    Teuchos::RCP<CompositeVectorSpace> dep_fac =
        Teuchos::rcp(new CompositeVectorSpace());
    dep_fac->SetOwned(false);
    dep_fac->SetGhosted(my_fac->Ghosted());
    dep_fac->SetMesh(my_fac->Mesh());
    dep_fac->AddComponent("cell", AmanziMesh::CELL, 1);

    // Loop over my dependencies, ensuring they meet the requirements.
    for (KeySet::const_iterator key=dependencies_.begin();
         key!=dependencies_.end(); ++key) {
      Teuchos::RCP<CompositeVectorSpace> fac = S->Require<CompositeVector,CompositeVectorSpace>(*key, Tags::NEXT);
      fac->Update(*dep_fac);
    }

    // Recurse into the tree to propagate info to leaves.
    for (KeySet::const_iterator key=dependencies_.begin();
         key!=dependencies_.end(); ++key) {
      S->RequireEvaluator(*key)->EnsureCompatibility(S);
    }
  }

}


// ---------------------------------------------------------------------------
// Updates the total derivative d(my_field)/d(wrt_field) in state S.
//
// This is derived to make sure that the derivatives created within this
// method only get cell quantities, not face quantites, in the same way that
// EnsureCompatibility() did.
// ---------------------------------------------------------------------------
void HeightEvaluator::UpdateFieldDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key) {

  Key dmy_key = Keys::getDerivKey(my_key_,wrt_key);
 
  Teuchos::RCP<CompositeVector> dmy;
  if (S->HasField(dmy_key)) {
    // Get the field...
    dmy = S->GetPtrW<CompositeVector>(dmy_key, my_key_);
  } else {
    // or create the field.  Note we have to do extra work that is normally
    // done by State in initialize.
    Teuchos::RCP<CompositeVectorSpace> my_fac = S->Require<CompositeVector,CompositeVectorSpace>(my_key_, Tags::NEXT);
    Teuchos::RCP<CompositeVectorSpace> new_fac =
      S->Require<CompositeVector,CompositeVectorSpace>(dmy_key, Tags::NEXT,  my_key_);
    new_fac->Update(*my_fac);
    dmy = Teuchos::rcp(new CompositeVector(*my_fac));
    S->SetData(dmy_key, my_key_, dmy);
    S->GetField(dmy_key,my_key_)->set_initialized();
    S->GetField(dmy_key,my_key_)->set_io_vis(false);
    S->GetField(dmy_key,my_key_)->set_io_checkpoint(false);
  }

  // Now update the values.
  dmy->PutScalar(0.0);

  // Only get derivatives in cells
  for (KeySet::const_iterator dep=dependencies_.begin();
       dep!=dependencies_.end(); ++dep) {
    Teuchos::RCP<CompositeVector> tmp = Teuchos::rcp(new CompositeVector(*dmy));

    if (wrt_key == *dep) {
      // partial F / partial x
      EvaluateFieldPartialDerivative_(S, wrt_key, tmp.ptr());
      dmy->ViewComponent("cell",false)->Update(1.0,
              *tmp->ViewComponent("cell",false), 1.0);

      if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
        *vo_->os() << dmy_key << " = " << (*tmp)("cell",0) << std::endl;
      }

    } else if (S->GetEvaluator(*dep)->IsDependency(S, wrt_key)) {
      // partial F / partial dep * ddep/dx
      // -- ddep/dx

      Key ddep_key = Keys::getDerivKey(*dep,wrt_key);

      Teuchos::RCP<const CompositeVector> ddep = S->GetPtr<CompositeVector>(ddep_key);
      // -- partial F / partial dep
      EvaluateFieldPartialDerivative_(S, *dep, tmp.ptr());

      if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
        *vo_->os() << ddep_key << " = " << (*ddep)("cell",0) << ", ";
        *vo_->os() << Keys::getDerivKey(my_key_, *dep) << " = " << (*tmp)("cell",0) << std::endl;
      }

      dmy->ViewComponent("cell",false)
          ->Multiply(1.0, *ddep->ViewComponent("cell",false),
                     *tmp->ViewComponent("cell",false), 1.0);
    }
  }

  if (dmy->HasComponent("face"))
    dmy->ViewComponent("face",false)->PutScalar(1.0);
}


void HeightEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
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
  const Epetra_MultiVector& rho = *S->GetPtr<CompositeVector>(dens_key_)
      ->ViewComponent("cell",false);

  const double& p_atm = *S->GetScalarData(patm_key_);
  const Epetra_Vector& gravity = *S->GetConstantVectorData(gravity_key_);
  double gz = -gravity[2];  // check this

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


void HeightEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {

  // -- cells need the function eval
  const Epetra_MultiVector& res_c = *result->ViewComponent("cell",false);
  const Epetra_MultiVector& pres_c = *S->GetPtr<CompositeVector>(pres_key_)
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& rho = *S->GetPtr<CompositeVector>(dens_key_)
      ->ViewComponent("cell",false);

  const double& p_atm = *S->GetScalarData(patm_key_);
  const Epetra_Vector& gravity = *S->GetConstantVectorData(gravity_key_);
  double gz = -gravity[2];  // check this

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
