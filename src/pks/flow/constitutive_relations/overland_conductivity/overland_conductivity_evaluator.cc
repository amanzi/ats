/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluates the conductivity of surface flow.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "Mesh_Algorithms.hh"
#include "independent_variable_field_evaluator.hh"
#include "overland_conductivity_evaluator.hh"
#include "manning_conductivity_model.hh"
#include "split_denominator_conductivity_model.hh"
#include "ponded_depth_passthrough_conductivity_model.hh"

namespace Amanzi {
namespace Flow {

OverlandConductivityEvaluator::OverlandConductivityEvaluator(Teuchos::ParameterList& plist)
    : SecondaryVariableFieldEvaluator(plist)
{
  Key domain = Keys::getDomain(my_key_);

  if (plist_.isParameter("height key") || plist_.isParameter("ponded depth key")
      || plist_.isParameter("height key suffix") || plist_.isParameter("ponded depth key suffix")) {
    Errors::Message message("OverlandConductivity: only use \"depth key\" or \"depth key suffix\", not \"height key\" or \"ponded depth key\".");
    Exceptions::amanzi_throw(message);
  }
  depth_key_ = Keys::readKey(plist_, domain, "depth", "ponded_depth");
  dependencies_.insert(depth_key_);

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


void
OverlandConductivityEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S)
{
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

    // Create a second factory with only cell entries
    CompositeVectorSpace cell_only_fac;
    cell_only_fac.SetMesh(my_fac->Mesh())
      ->SetGhosted()->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

    // Loop over my dependencies, ensuring they meet the requirements.
    for (const auto& key : dependencies_) {
      if (key == my_key_) {
        Errors::Message msg;
        msg << "Evaluator for key \"" << my_key_ << "\" depends upon itself.";
        Exceptions::amanzi_throw(msg);
      } else if (key == slope_key_ || key == coef_key_) {
        Teuchos::RCP<CompositeVectorSpace> fac = S->RequireField(key);
        fac->Update(cell_only_fac);
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


// Required methods from SecondaryVariableFieldEvaluator
void OverlandConductivityEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
  Teuchos::RCP<const CompositeVector> depth = S->GetFieldData(depth_key_);
  Teuchos::RCP<const CompositeVector> slope = S->GetFieldData(slope_key_);
  Teuchos::RCP<const CompositeVector> coef = S->GetFieldData(coef_key_);

  { // cell component
    const Epetra_MultiVector& depth_v = *depth->ViewComponent("cell",false);
    const Epetra_MultiVector& slope_v = *slope->ViewComponent("cell",false);
    const Epetra_MultiVector& coef_v = *coef->ViewComponent("cell",false);

#ifdef ENABLE_DBC
  double min_coef = 1.;
  coef_v.MinValue(&min_coef);
  if (min_coef <= 0.) {
    Errors::Message message("Overland Conductivity Evaluator: Manning coeficient has at least one value that is non-positive.  Regions do not cover or do not set component: \"cell\"");
    Exceptions::amanzi_throw(message);
  }
#endif

    Epetra_MultiVector& result_v = *result->ViewComponent("cell",false);

    int ncomp = result->size("cell", false);
    if (dt_swe_factor_ > 0) {
      for (int i=0; i!=ncomp; ++i) {
        double new_snow = dt_swe_factor_ * depth_v[0][i];
        result_v[0][i] = model_->Conductivity(new_snow, slope_v[0][i], coef_v[0][i]);
      }
    } else {
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->Conductivity(depth_v[0][i], slope_v[0][i], coef_v[0][i]);
      }
    }

    if (dens_) {
      const Epetra_MultiVector& dens_v = *S->GetFieldData(dens_key_)->ViewComponent("cell",false);
      for (int i=0; i!=ncomp; ++i) result_v[0][i] *= dens_v[0][i];
    }
  }

  if (result->HasComponent("boundary_face")) {
    const Epetra_MultiVector& depth_v = *depth->ViewComponent("boundary_face",false);
    const Epetra_MultiVector& slope_v = *slope->ViewComponent("cell",false);
    const Epetra_MultiVector& coef_v = *coef->ViewComponent("cell",false);
    Epetra_MultiVector& result_v = *result->ViewComponent("boundary_face",false);
    const auto& mesh = *(result->Mesh());

    int ncomp = result->size("boundary_face", false);
    if (dt_swe_factor_ > 0) {
      for (int i=0; i!=ncomp; ++i) {
        AmanziMesh::Entity_ID c = AmanziMesh::getBoundaryFaceInternalCell(mesh, i);
        double new_snow = dt_swe_factor_ * depth_v[0][i];
        result_v[0][i] = model_->Conductivity(new_snow, slope_v[0][c], coef_v[0][c]);
      }
    } else {
      for (int i=0; i!=ncomp; ++i) {
        AmanziMesh::Entity_ID c = AmanziMesh::getBoundaryFaceInternalCell(mesh, i);
        result_v[0][i] = model_->Conductivity(depth_v[0][i], slope_v[0][c], coef_v[0][c]);
      }
    }

    if (dens_) {
      const Epetra_MultiVector& dens_v = *S->GetFieldData(dens_key_)->ViewComponent("boundary_face",false);
      for (int i=0; i!=ncomp; ++i) result_v[0][i] *= dens_v[0][i];
    }
  }
}


void
OverlandConductivityEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  Teuchos::RCP<const CompositeVector> depth = S->GetFieldData(depth_key_);
  Teuchos::RCP<const CompositeVector> slope = S->GetFieldData(slope_key_);
  Teuchos::RCP<const CompositeVector> coef = S->GetFieldData(coef_key_);

  if (wrt_key == depth_key_) {
    { // cell component
      const Epetra_MultiVector& depth_v = *depth->ViewComponent("cell",false);
      const Epetra_MultiVector& slope_v = *slope->ViewComponent("cell",false);
      const Epetra_MultiVector& coef_v = *coef->ViewComponent("cell",false);
      Epetra_MultiVector& result_v = *result->ViewComponent("cell",false);

      int ncomp = result->size("cell", false);
      if (dt_swe_factor_ > 0.) {
        for (int i=0; i!=ncomp; ++i) {
          double new_snow = dt_swe_factor_ * depth_v[0][i];
          result_v[0][i] = model_->DConductivityDDepth(new_snow, slope_v[0][i], coef_v[0][i])
                           * dt_swe_factor_;
        }
      } else {
        for (int i=0; i!=ncomp; ++i) {
          result_v[0][i] = model_->DConductivityDDepth(depth_v[0][i], slope_v[0][i], coef_v[0][i]);
        }
      }

      if (dens_) {
        const Epetra_MultiVector& dens_v = *S->GetFieldData(dens_key_)->ViewComponent("cell",false);
        for (int i=0; i!=ncomp; ++i) {
          result_v[0][i] *= dens_v[0][i];
        }
      }
    }

    if (result->HasComponent("boundary_face")) {
      const Epetra_MultiVector& depth_v = *depth->ViewComponent("boundary_face",false);
      const Epetra_MultiVector& slope_v = *slope->ViewComponent("cell",false);
      const Epetra_MultiVector& coef_v = *coef->ViewComponent("cell",false);
      Epetra_MultiVector& result_v = *result->ViewComponent("boundary_face",false);
      const auto& mesh = *result->Mesh();

      int ncomp = result->size("boundary_face", false);
      if (dt_swe_factor_ > 0) {
        for (int i=0; i!=ncomp; ++i) {
          AmanziMesh::Entity_ID c = AmanziMesh::getBoundaryFaceInternalCell(mesh, i);
          double new_snow = dt_swe_factor_ * depth_v[0][i];
          result_v[0][i] = model_->DConductivityDDepth(new_snow, slope_v[0][c], coef_v[0][c])
            * dt_swe_factor_;
        }
      } else {
        for (int i=0; i!=ncomp; ++i) {
          AmanziMesh::Entity_ID c = AmanziMesh::getBoundaryFaceInternalCell(mesh, i);
          result_v[0][i] = model_->DConductivityDDepth(depth_v[0][i], slope_v[0][c], coef_v[0][c]);
        }
      }

      if (dens_) {
        const Epetra_MultiVector& dens_v = *S->GetFieldData(dens_key_)->ViewComponent("boundary_face",false);
        for (int i=0; i!=ncomp; ++i) result_v[0][i] *= dens_v[0][i];
      }
    }

  } else if (wrt_key == dens_key_) {
    AMANZI_ASSERT(dens_);
    { // cell component
      const Epetra_MultiVector& depth_v = *depth->ViewComponent("cell",false);
      const Epetra_MultiVector& slope_v = *slope->ViewComponent("cell",false);
      const Epetra_MultiVector& coef_v = *coef->ViewComponent("cell",false);
      Epetra_MultiVector& result_v = *result->ViewComponent("cell",false);

      int ncomp = result->size("cell", false);
      if (dt_swe_factor_ > 0.) {
        for (int i=0; i!=ncomp; ++i) {
          double new_snow = dt_swe_factor_ * depth_v[0][i];
          result_v[0][i] = model_->Conductivity(new_snow, slope_v[0][i], coef_v[0][i]);
        }
      } else {
        for (int i=0; i!=ncomp; ++i) {
          result_v[0][i] = model_->Conductivity(depth_v[0][i], slope_v[0][i], coef_v[0][i]);
        }
      }
    }

    if (result->HasComponent("boundary_face")) {
      const Epetra_MultiVector& depth_v = *depth->ViewComponent("boundary_face",false);
      const Epetra_MultiVector& slope_v = *slope->ViewComponent("cell",false);
      const Epetra_MultiVector& coef_v = *coef->ViewComponent("cell",false);
      Epetra_MultiVector& result_v = *result->ViewComponent("boundary_face",false);
      const auto& mesh = *result->Mesh();

      int ncomp = result->size("boundary_face", false);
      if (dt_swe_factor_ > 0) {
        for (int i=0; i!=ncomp; ++i) {
          AmanziMesh::Entity_ID c = AmanziMesh::getBoundaryFaceInternalCell(mesh, i);
          double new_snow = dt_swe_factor_ * depth_v[0][i];
          result_v[0][i] = model_->Conductivity(new_snow, slope_v[0][c], coef_v[0][c]);
        }
      } else {
        for (int i=0; i!=ncomp; ++i) {
          AmanziMesh::Entity_ID c = AmanziMesh::getBoundaryFaceInternalCell(mesh, i);
          result_v[0][i] = model_->Conductivity(depth_v[0][i], slope_v[0][c], coef_v[0][c]);
        }
      }
    }
  } else {
    // FIX ME -- need to add derivatives of conductivity model wrt slope, coef --etc
    result->PutScalar(0.);
  }
}


} //namespace
} //namespace

