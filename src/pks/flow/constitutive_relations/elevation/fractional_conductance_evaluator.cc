/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#include "Mesh_Algorithms.hh"
#include "fractional_conductance_evaluator.hh"
#include "subgrid_microtopography.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

FractionalConductanceEvaluator::FractionalConductanceEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  Key domain = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  vpd_key_ = Keys::readKey(plist_, domain, "volumetric ponded depth", "volumetric_ponded_depth");
  dependencies_.insert(KeyTag{vpd_key_, tag});

  mobile_depth_key_ = Keys::readKey(plist_, domain, "mobile depth", "mobile_depth");
  dependencies_.insert(KeyTag{mobile_depth_key_, tag});

  delta_max_key_ = Keys::readKey(plist_, domain, "microtopographic relief", "microtopographic_relief");
  dependencies_.insert(KeyTag{delta_max_key_, tag});

  delta_ex_key_ = Keys::readKey(plist_, domain, "excluded volume", "excluded_volume");
  dependencies_.insert(KeyTag{delta_ex_key_, tag});

  depr_depth_key_ = Keys::readKey(plist_, domain, "depression depth", "depression_depth");
  dependencies_.insert(KeyTag{depr_depth_key_, tag});
}


Teuchos::RCP<Evaluator>
FractionalConductanceEvaluator::Clone() const {
  return Teuchos::rcp(new FractionalConductanceEvaluator(*this));
}


void FractionalConductanceEvaluator::Evaluate_(const State& S,
        const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> vpd_v = S->GetPtr<CompositeVector>(vpd_key_, tag);
  Teuchos::RCP<const CompositeVector> del_max_v = S->GetPtr<CompositeVector>(delta_max_key_, tag);
  Teuchos::RCP<const CompositeVector> del_ex_v = S->GetPtr<CompositeVector>(delta_ex_key_, tag);
  Teuchos::RCP<const CompositeVector> depr_depth_v = S->GetPtr<CompositeVector>(depr_depth_key_, tag);
  Teuchos::RCP<const CompositeVector> mobile_depth_v = S->GetPtr<CompositeVector>(mobile_depth_key_, tag);
  const auto& mesh = *result->Mesh();

  for (const auto& comp : *result) {
    AMANZI_ASSERT(comp == "cell" || comp == "boundary_face");
    bool is_internal_comp = comp == "boundary_face";
    Key internal_comp = is_internal_comp ? "cell" : comp;

    const auto& vpd = *vpd_v->ViewComponent(comp,false);
    const auto& mobile_depth = *mobile_depth_v->ViewComponent(comp,false);

    const auto& del_max = *del_max_v->ViewComponent(internal_comp,false);
    const auto& del_ex = *del_ex_v->ViewComponent(internal_comp,false);
    const auto& depr_depth = *depr_depth_v->ViewComponent(internal_comp,false);
    auto& res = *result->ViewComponent(comp,false);

    int ncomp = result->size(comp, false);
    for (int i=0; i!=ncomp; ++i) {
      int ii = is_internal_comp ? AmanziMesh::getBoundaryFaceInternalCell(mesh, i) : i;
      if (mobile_depth[0][i] <= 0.0) {
        res[0][i] = 0;
      } else {
        double vol_depr_depth = Microtopography::volumetricDepth(depr_depth[0][ii],
                del_max[0][ii], del_ex[0][ii]);
        res[0][i] = (vpd[0][i] - vol_depr_depth) / (mobile_depth[0][i]);
      }
    }
  }
}


void
FractionalConductanceEvaluator::EvaluatePartialDerivative_(const State& S,
        const Key& wrt_key, const Tag& wrt_tag, const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> vpd_v = S->GetPtr<CompositeVector>(vpd_key_, tag);
  Teuchos::RCP<const CompositeVector> del_max_v = S->GetPtr<CompositeVector>(delta_max_key_, tag);
  Teuchos::RCP<const CompositeVector> del_ex_v = S->GetPtr<CompositeVector>(delta_ex_key_, tag);
  Teuchos::RCP<const CompositeVector> depr_depth_v = S->GetPtr<CompositeVector>(depr_depth_key_, tag);
  Teuchos::RCP<const CompositeVector> mobile_depth_v = S->GetPtr<CompositeVector>(mobile_depth_key_, tag);
  const auto& mesh = *result->Mesh();

  if (wrt_key == mobile_depth_key_) {
    for (const auto& comp : *result) {
      AMANZI_ASSERT(comp == "cell" || comp == "boundary_face");
      bool is_internal_comp = comp == "boundary_face";
      Key internal_comp = is_internal_comp ? "cell" : comp;

      const auto& vpd = *vpd_v->ViewComponent(comp,false);
      const auto& mobile_depth = *mobile_depth_v->ViewComponent(comp,false);

      const auto& del_max = *del_max_v->ViewComponent(internal_comp,false);
      const auto& del_ex = *del_ex_v->ViewComponent(internal_comp,false);
      const auto& depr_depth = *depr_depth_v->ViewComponent(internal_comp,false);
      auto& res = *result->ViewComponent(comp,false);

      int ncomp = result->size(comp, false);
      for (int i=0; i!=ncomp; ++i) {
        int ii = is_internal_comp ? AmanziMesh::getBoundaryFaceInternalCell(mesh, i) : i;
        if (mobile_depth[0][i] <= 0.0) {
          res[0][i] = 0;
        } else {
          double vol_depr_depth = Microtopography::volumetricDepth(depr_depth[0][ii], del_max[0][ii], del_ex[0][ii]);
          res[0][i] = - (vpd[0][i] - vol_depr_depth) * std::pow(mobile_depth[0][ii],-2.);
        }
      }
    }
  } else if (wrt_key == vpd_key_) {
    for (const auto& comp : *result) {
      const auto& mobile_depth = *mobile_depth_v->ViewComponent(comp,false);
      auto& res = *result->ViewComponent(comp,false);

      int ncomp = result->size(comp, false);
      for (int i=0; i!=ncomp; ++i) {
        if (mobile_depth[0][i] <= 0.0) {
          res[0][i] = 0;
        } else {
          res[0][i] = 1./mobile_depth[0][i];
        }
      }
    }
  } else {
    Errors::Message msg("FractionalConductanceEvaluator: Not Implemented: no derivatives implemented other than mobile depth and volumentric ponded depth.");
    Exceptions::amanzi_throw(msg);
  }
}


void FractionalConductanceEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S)
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

      if (key == depr_depth_key_ ||
          key == delta_ex_key_ ||
          key == delta_max_key_) {
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
} //namespace
