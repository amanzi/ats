/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluates the conductivity of surface flow.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "Mesh_Algorithms.hh"
#include "overland_conductivity_evaluator.hh"
#include "manning_conductivity_model.hh"

namespace Amanzi {
namespace Flow {

OverlandConductivityEvaluator::OverlandConductivityEvaluator(Teuchos::ParameterList& plist)
    : SecondaryVariableFieldEvaluator(plist)
{
  Key domain = Keys::getDomain(my_key_);

  if (plist_.isParameter("height key") ||
      plist_.isParameter("ponded depth key") ||
      plist_.isParameter("depth key") ||
      plist_.isParameter("height key suffix") ||
      plist_.isParameter("ponded depth key suffix") ||
      plist_.isParameter("depth key suffix")) {
    Errors::Message message("OverlandConductivity: only use \"mobile depth key\" or \"mobile depth key suffix\", not \"height key\" or \"ponded depth key\" or \"depth key\".");
    Exceptions::amanzi_throw(message);
  }
  mobile_depth_key_ = Keys::readKey(plist_, domain, "mobile depth", "ponded_depth");
  dependencies_.insert(mobile_depth_key_);

  slope_key_ = Keys::readKey(plist_, domain, "slope", "slope_magnitude");
  dependencies_.insert(slope_key_);

  coef_key_ = Keys::readKey(plist_, domain, "coefficient", "manning_coefficient");
  dependencies_.insert(coef_key_);

  dt_swe_factor_ = plist_.get<double>("dt factor [s]", -1);
  if (dt_swe_factor_ > 0) {
    double swe_factor = plist_.get<double>("swe density factor [-]", 10.0);
    dt_swe_factor_ *= swe_factor;
  }

  dens_ = plist_.get<bool>("include density", true);
  if (dens_) {
    dens_key_ = Keys::readKey(plist_, domain, "molar density liquid", "molar_density_liquid");
    dependencies_.insert(dens_key_);
  }

  // create the model
  Teuchos::ParameterList& sublist = plist_.sublist("overland conductivity model");
  model_ = Teuchos::rcp(new ManningConductivityModel(sublist));
}


Teuchos::RCP<FieldEvaluator>
OverlandConductivityEvaluator::Clone() const
{
  return Teuchos::rcp(new OverlandConductivityEvaluator(*this));
}


// Required methods from SecondaryVariableFieldEvaluator
void OverlandConductivityEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
  Teuchos::RCP<const CompositeVector> depth = S->GetFieldData(mobile_depth_key_);
  Teuchos::RCP<const CompositeVector> slope = S->GetFieldData(slope_key_);
  Teuchos::RCP<const CompositeVector> coef = S->GetFieldData(coef_key_);

#ifdef ENABLE_DBC
  double min_coef = 1.;
  coef->MinValue(&min_coef);
  if (min_coef <= 1.e-12) {
    Errors::Message message("Overland Conductivity Evaluator: Manning coeficient has at least one value that is non-positive.  Perhaps you forgot to set the \"boundary_face\" component?");
    Exceptions::amanzi_throw(message);
  }
#endif

  // Note, the logic here splits vectors into two groups -- those that will be
  // defined on all components and those that do not make sense to define on
  // not-cells.  These vectors are used in a wierd way on boundary faces,
  // getting the internal cell and using that.  This currently cannot be used
  // on any other components than "cell" and "boundary_face".
  const AmanziMesh::Mesh& mesh = *result->Mesh();
  for (const auto& comp : *result) {
    AMANZI_ASSERT(comp == "cell" || comp == "boundary_face");
    bool is_internal_comp = comp == "boundary_face";
    Key internal_comp = is_internal_comp ? "cell" : comp;

    const Epetra_MultiVector& slope_v = *slope->ViewComponent(internal_comp,false);
    const Epetra_MultiVector& coef_v = *coef->ViewComponent(internal_comp,false);

    const Epetra_MultiVector& depth_v = *depth->ViewComponent(comp,false);
    Epetra_MultiVector& result_v = *result->ViewComponent(comp,false);

    int ncomp = result->size(comp, false);
    if (dt_swe_factor_ > 0) {
      for (int i=0; i!=ncomp; ++i) {
        int ii = is_internal_comp ? AmanziMesh::getBoundaryFaceInternalCell(mesh, i) : i;
        double new_snow = dt_swe_factor_ * depth_v[0][i];
        result_v[0][i] = model_->Conductivity(new_snow, slope_v[0][ii], coef_v[0][ii]);
      }
    } else {
      for (int i=0; i!=ncomp; ++i) {
        int ii = is_internal_comp ? AmanziMesh::getBoundaryFaceInternalCell(mesh, i) : i;
        result_v[0][i] = model_->Conductivity(depth_v[0][i], slope_v[0][ii], coef_v[0][ii]);
      }
    }

    if (dens_) {
      const Epetra_MultiVector& dens_v = *S->GetFieldData(dens_key_)->ViewComponent(comp,false);
      for (int i=0; i!=ncomp; ++i) result_v[0][i] *= dens_v[0][i];
    }
  }
}


void
OverlandConductivityEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  // NOTE, we can only differentiate with respect to quantities that exist on
  // all entities, not just cell entities.
  Teuchos::RCP<const CompositeVector> depth = S->GetFieldData(mobile_depth_key_);
  Teuchos::RCP<const CompositeVector> slope = S->GetFieldData(slope_key_);
  Teuchos::RCP<const CompositeVector> coef = S->GetFieldData(coef_key_);
  const AmanziMesh::Mesh& mesh = *result->Mesh();

  if (wrt_key == mobile_depth_key_) {
    for (const auto& comp : *result) {
      AMANZI_ASSERT(comp == "cell" || comp == "boundary_face");
      bool is_internal_comp = comp == "boundary_face";
      Key internal_comp = is_internal_comp ? "cell" : comp;

      const Epetra_MultiVector& slope_v = *slope->ViewComponent(internal_comp,false);
      const Epetra_MultiVector& coef_v = *coef->ViewComponent(internal_comp,false);

      const Epetra_MultiVector& depth_v = *depth->ViewComponent(comp,false);
      Epetra_MultiVector& result_v = *result->ViewComponent(comp,false);

      int ncomp = result->size(comp, false);
      if (dt_swe_factor_ > 0.) {
        for (int i=0; i!=ncomp; ++i) {
          int ii = is_internal_comp ? AmanziMesh::getBoundaryFaceInternalCell(mesh, i) : i;
          double new_snow = dt_swe_factor_ * depth_v[0][i];
          result_v[0][i] = model_->DConductivityDDepth(new_snow, slope_v[0][ii], coef_v[0][ii])
                           * dt_swe_factor_;
        }
      } else {
        for (int i=0; i!=ncomp; ++i) {
          int ii = is_internal_comp ? AmanziMesh::getBoundaryFaceInternalCell(mesh, i) : i;
          result_v[0][i] = model_->DConductivityDDepth(depth_v[0][i], slope_v[0][ii], coef_v[0][ii]);
        }
      }

      if (dens_) {
        const Epetra_MultiVector& dens_v = *S->GetFieldData(dens_key_)->ViewComponent(comp,false);
        for (int i=0; i!=ncomp; ++i) {
          result_v[0][i] *= dens_v[0][i];
        }
      }
    }

  } else if (wrt_key == dens_key_) {
    AMANZI_ASSERT(dens_);
    for (const auto& comp : *result) {
      AMANZI_ASSERT(comp == "cell" || comp == "boundary_face");
      bool is_internal_comp = comp == "boundary_face";
      Key internal_comp = is_internal_comp ? "cell" : comp;

      const Epetra_MultiVector& slope_v = *slope->ViewComponent(internal_comp,false);
      const Epetra_MultiVector& coef_v = *coef->ViewComponent(internal_comp,false);

      const Epetra_MultiVector& depth_v = *depth->ViewComponent(comp,false);
      Epetra_MultiVector& result_v = *result->ViewComponent(comp,false);

      int ncomp = result->size(comp, false);
      if (dt_swe_factor_ > 0.) {
        for (int i=0; i!=ncomp; ++i) {
          int ii = is_internal_comp ? AmanziMesh::getBoundaryFaceInternalCell(mesh, i) : i;
          double new_snow = dt_swe_factor_ * depth_v[0][i];
          result_v[0][i] = model_->Conductivity(new_snow, slope_v[0][ii], coef_v[0][ii]);
        }
      } else {
        for (int i=0; i!=ncomp; ++i) {
          int ii = is_internal_comp ? AmanziMesh::getBoundaryFaceInternalCell(mesh, i) : i;
          result_v[0][i] = model_->Conductivity(depth_v[0][i], slope_v[0][ii], coef_v[0][ii]);
        }
      }
    }

  } else {
    result->PutScalar(0.);
  }
}


void OverlandConductivityEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S)
{
  // Ensure my field exists.  Requirements should be already set.
  AMANZI_ASSERT(my_key_ != std::string(""));
  Teuchos::RCP<CompositeVectorSpace> my_fac = S->RequireField(my_key_, my_key_);

  // check plist for vis or checkpointing control
  bool io_my_key = plist_.get<bool>("visualize", true);
  S->GetField(my_key_, my_key_)->set_io_vis(io_my_key);
  bool checkpoint_my_key = plist_.get<bool>("checkpoint", false);
  S->GetField(my_key_, my_key_)->set_io_checkpoint(checkpoint_my_key);

  // If my requirements have not yet been set, we'll have to hope they
  // get set by someone later.  For now just defer.
  if (my_fac->Mesh() != Teuchos::null) {
    // Create an unowned factory to check my dependencies.
    Teuchos::RCP<CompositeVectorSpace> dep_fac =
        Teuchos::rcp(new CompositeVectorSpace(*my_fac));
    dep_fac->SetOwned(false);

    Teuchos::RCP<CompositeVectorSpace> no_bf_dep_fac;
    if (dep_fac->HasComponent("boundary_face")) {
      no_bf_dep_fac = Teuchos::rcp(new CompositeVectorSpace());
      no_bf_dep_fac->SetMesh(dep_fac->Mesh())
        ->SetGhosted(true)
        ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    } else {
      no_bf_dep_fac = dep_fac;
    }

    // Loop over my dependencies, ensuring they meet the requirements.
    for (const auto& key : dependencies_) {
      if (key == my_key_) {
        Errors::Message msg;
        msg << "Evaluator for key \"" << my_key_ << "\" depends upon itself.";
        Exceptions::amanzi_throw(msg);
      }

      if (key == coef_key_ || key == slope_key_) {
        Teuchos::RCP<CompositeVectorSpace> fac = S->RequireField(key);
        fac->Update(*no_bf_dep_fac);
      } else {
        Teuchos::RCP<CompositeVectorSpace> fac = S->RequireField(key);
        fac->Update(*dep_fac);
      }
    }

    // Recurse into the tree to propagate info to leaves.
    for (const auto& key : dependencies_) {
      S->RequireFieldEvaluator(key)->EnsureCompatibility(S);
    }
  }
}



} //namespace
} //namespace

