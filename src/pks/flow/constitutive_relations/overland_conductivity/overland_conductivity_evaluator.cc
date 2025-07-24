/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Evaluates the conductivity of surface flow.

*/

#include "MeshAlgorithms.hh"
#include "overland_conductivity_evaluator.hh"
#include "manning_conductivity_model.hh"

namespace Amanzi {
namespace Flow {

OverlandConductivityEvaluator::OverlandConductivityEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  Key domain = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  if (plist_.isParameter("height key") || plist_.isParameter("ponded depth key") ||
      plist_.isParameter("depth key") || plist_.isParameter("height key suffix") ||
      plist_.isParameter("ponded depth key suffix") || plist_.isParameter("depth key suffix")) {
    Errors::Message message(
      "OverlandConductivity: only use \"mobile depth key\" or \"mobile depth key suffix\", not "
      "\"height key\" or \"ponded depth key\" or \"depth key\".");
    Exceptions::amanzi_throw(message);
  }


  mobile_depth_key_ = Keys::readKey(plist_, domain, "mobile depth", "ponded_depth");
  dependencies_.insert(KeyTag{ mobile_depth_key_, tag });

  slope_key_ = Keys::readKey(plist_, domain, "slope", "slope_magnitude");
  dependencies_.insert(KeyTag{ slope_key_, tag });

  coef_key_ = Keys::readKey(plist_, domain, "coefficient", "manning_coefficient");
  dependencies_.insert(KeyTag{ coef_key_, tag });

  dt_swe_factor_ = plist_.get<double>("dt factor [s]", -1);
  if (dt_swe_factor_ > 0) {
    double swe_factor = plist_.get<double>("swe density factor [-]", 10.0);
    dt_swe_factor_ *= swe_factor;
  }

  dens_ = plist_.get<bool>("include density", true);
  if (dens_) {
    dens_key_ = Keys::readKey(plist_, domain, "molar density liquid", "molar_density_liquid");
    dependencies_.insert(KeyTag{ dens_key_, tag });
  }

  // create the model
  Teuchos::ParameterList& sublist = plist_.sublist("overland conductivity model");
  model_ = Teuchos::rcp(new ManningConductivityModel(sublist));
}


Teuchos::RCP<Evaluator>
OverlandConductivityEvaluator::Clone() const
{
  return Teuchos::rcp(new OverlandConductivityEvaluator(*this));
}


// Required methods from EvaluatorSecondaryMonotypeCV
void
OverlandConductivityEvaluator::Evaluate_(const State& S,
                                         const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> depth = S.GetPtr<CompositeVector>(mobile_depth_key_, tag);
  Teuchos::RCP<const CompositeVector> slope = S.GetPtr<CompositeVector>(slope_key_, tag);
  Teuchos::RCP<const CompositeVector> coef = S.GetPtr<CompositeVector>(coef_key_, tag);

#ifdef ENABLE_DBC
  double min_coef = 1.;
  coef->MinValue(&min_coef);
  if (min_coef <= 1.e-12) {
    Errors::Message message(
      "Overland Conductivity Evaluator: Manning coeficient has at least one value that is "
      "non-positive.  Perhaps you forgot to set the \"boundary_face\" component?");
    Exceptions::amanzi_throw(message);
  }
#endif

  // Note, the logic here splits vectors into two groups -- those that will be
  // defined on all components and those that do not make sense to define on
  // not-cells.  These vectors are used in a wierd way on boundary faces,
  // getting the internal cell and using that.  This currently cannot be used
  // on any other components than "cell" and "boundary_face".
  const AmanziMesh::Mesh& mesh = *result[0]->Mesh();
  for (const auto& comp : *result[0]) {
    AMANZI_ASSERT(comp == "cell" || comp == "boundary_face");
    bool is_internal_comp = comp == "boundary_face";
    Key internal_comp = is_internal_comp ? "cell" : comp;

    const Epetra_MultiVector& slope_v = *slope->ViewComponent(internal_comp, false);
    const Epetra_MultiVector& coef_v = *coef->ViewComponent(internal_comp, false);

    const Epetra_MultiVector& depth_v = *depth->ViewComponent(comp, false);
    Epetra_MultiVector& result_v = *result[0]->ViewComponent(comp, false);

    int ncomp = result[0]->size(comp, false);
    if (dt_swe_factor_ > 0) {
      for (int i = 0; i != ncomp; ++i) {
        int ii = is_internal_comp ? AmanziMesh::getBoundaryFaceInternalCell(mesh, i) : i;
        double new_snow = dt_swe_factor_ * depth_v[0][i];
        result_v[0][i] = model_->Conductivity(new_snow, slope_v[0][ii], coef_v[0][ii]);
      }
    } else {
      for (int i = 0; i != ncomp; ++i) {
        int ii = is_internal_comp ? AmanziMesh::getBoundaryFaceInternalCell(mesh, i) : i;
        result_v[0][i] = model_->Conductivity(depth_v[0][i], slope_v[0][ii], coef_v[0][ii]);
      }
    }

    if (dens_) {
      const Epetra_MultiVector& dens_v =
        *S.Get<CompositeVector>(dens_key_, tag).ViewComponent(comp, false);
      for (int i = 0; i != ncomp; ++i) result_v[0][i] *= dens_v[0][i];
    }
  }
}


void
OverlandConductivityEvaluator::EvaluatePartialDerivative_(
  const State& S,
  const Key& wrt_key,
  const Tag& wrt_tag,
  const std::vector<CompositeVector*>& result)
{
  // NOTE, we can only differentiate with respect to quantities that exist on
  // all entities, not just cell entities.
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> depth = S.GetPtr<CompositeVector>(mobile_depth_key_, tag);
  Teuchos::RCP<const CompositeVector> slope = S.GetPtr<CompositeVector>(slope_key_, tag);
  Teuchos::RCP<const CompositeVector> coef = S.GetPtr<CompositeVector>(coef_key_, tag);

  const AmanziMesh::Mesh& mesh = *result[0]->Mesh();

  if (wrt_key == mobile_depth_key_) {
    for (const auto& comp : *result[0]) {
      AMANZI_ASSERT(comp == "cell" || comp == "boundary_face");
      bool is_internal_comp = comp == "boundary_face";
      Key internal_comp = is_internal_comp ? "cell" : comp;

      const Epetra_MultiVector& slope_v = *slope->ViewComponent(internal_comp, false);
      const Epetra_MultiVector& coef_v = *coef->ViewComponent(internal_comp, false);

      const Epetra_MultiVector& depth_v = *depth->ViewComponent(comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(comp, false);

      int ncomp = result[0]->size(comp, false);
      if (dt_swe_factor_ > 0.) {
        for (int i = 0; i != ncomp; ++i) {
          int ii = is_internal_comp ? AmanziMesh::getBoundaryFaceInternalCell(mesh, i) : i;
          double new_snow = dt_swe_factor_ * depth_v[0][i];
          result_v[0][i] =
            model_->DConductivityDDepth(new_snow, slope_v[0][ii], coef_v[0][ii]) * dt_swe_factor_;
        }
      } else {
        for (int i = 0; i != ncomp; ++i) {
          int ii = is_internal_comp ? AmanziMesh::getBoundaryFaceInternalCell(mesh, i) : i;
          result_v[0][i] =
            model_->DConductivityDDepth(depth_v[0][i], slope_v[0][ii], coef_v[0][ii]);
        }
      }

      if (dens_) {
        const Epetra_MultiVector& dens_v =
          *S.Get<CompositeVector>(dens_key_, tag).ViewComponent(comp, false);
        for (int i = 0; i != ncomp; ++i) {
          result_v[0][i] *= dens_v[0][i];
        }
      }
    }

  } else if (wrt_key == dens_key_) {
    AMANZI_ASSERT(dens_);

    for (const auto& comp : *result[0]) {
      AMANZI_ASSERT(comp == "cell" || comp == "boundary_face");
      bool is_internal_comp = comp == "boundary_face";
      Key internal_comp = is_internal_comp ? "cell" : comp;

      const Epetra_MultiVector& slope_v = *slope->ViewComponent(internal_comp, false);
      const Epetra_MultiVector& coef_v = *coef->ViewComponent(internal_comp, false);

      const Epetra_MultiVector& depth_v = *depth->ViewComponent(comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(comp, false);

      int ncomp = result[0]->size(comp, false);
      if (dt_swe_factor_ > 0.) {
        for (int i = 0; i != ncomp; ++i) {
          int ii = is_internal_comp ? AmanziMesh::getBoundaryFaceInternalCell(mesh, i) : i;
          double new_snow = dt_swe_factor_ * depth_v[0][i];
          result_v[0][i] = model_->Conductivity(new_snow, slope_v[0][ii], coef_v[0][ii]);
        }
      } else {
        for (int i = 0; i != ncomp; ++i) {
          int ii = is_internal_comp ? AmanziMesh::getBoundaryFaceInternalCell(mesh, i) : i;
          result_v[0][i] = model_->Conductivity(depth_v[0][i], slope_v[0][ii], coef_v[0][ii]);
        }
      }
    }

  } else {
    // FIX ME -- need to add derivatives of conductivity model wrt slope, coef --etc
    result[0]->PutScalar(0.);
  }
}


void
OverlandConductivityEvaluator::EnsureCompatibility_ToDeps_(State& S)
{
  // Ensure my field exists.  Requirements should be already set.
  const auto& my_fac = S.Require<CompositeVector, CompositeVectorSpace>(my_keys_.front().first,
                                                                        my_keys_.front().second);

  // If my requirements have not yet been set, we'll have to hope they
  // get set by someone later.  For now just defer.
  if (my_fac.Mesh() != Teuchos::null) {
    // Create an unowned factory to check my dependencies.
    Teuchos::RCP<CompositeVectorSpace> dep_fac = Teuchos::rcp(new CompositeVectorSpace(my_fac));
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
    for (const auto& key_tag : dependencies_) {
      if (key_tag.first == coef_key_ || key_tag.first == slope_key_) {
        S.Require<CompositeVector, CompositeVectorSpace>(key_tag.first, key_tag.second)
          .Update(*no_bf_dep_fac);
      } else {
        S.Require<CompositeVector, CompositeVectorSpace>(key_tag.first, key_tag.second)
          .Update(*dep_fac);
      }
    }
  }
}


} // namespace Flow
} // namespace Amanzi
