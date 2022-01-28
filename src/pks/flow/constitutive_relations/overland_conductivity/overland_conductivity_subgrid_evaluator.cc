/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluates the conductivity of surface flow subgrid model.

  Authors: Ahmad Jan (jana@ornl.gov)
*/

#include "overland_conductivity_subgrid_evaluator.hh"
#include "manning_conductivity_model.hh"

namespace Amanzi {
namespace Flow {

OverlandConductivitySubgridEvaluator::OverlandConductivitySubgridEvaluator(Teuchos::ParameterList& plist)
    : EvaluatorSecondaryMonotypeCV(plist)
{
  Key domain = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  if (plist_.isParameter("height key") ||
      plist_.isParameter("ponded depth key") ||
      plist_.isParameter("depth key") ||
      plist_.isParameter("height key suffix") ||
      plist_.isParameter("ponded depth key suffix") ||
      plist_.isParameter("depth key suffix")) {
    Errors::Message message("OverlandConductivitySubgrid: only use \"mobile depth key\" or \"mobile depth key suffix\", not \"height key\" or \"ponded depth key\" or \"depth key\".");
    Exceptions::amanzi_throw(message);
  }

  mobile_depth_key_ = Keys::readKey(plist_, domain, "mobile depth", "mobile_depth");
  dependencies_.insert(KeyTag{mobile_depth_key_, tag});

  slope_key_ = Keys::readKey(plist_, domain, "slope", "slope_magnitude");
  dependencies_.insert(KeyTag{slope_key_, tag});

  coef_key_ = Keys::readKey(plist_, domain, "coefficient", "manning_coefficient");
  dependencies_.insert(KeyTag{coef_key_, tag});

  dens_key_ = Keys::readKey(plist_, domain, "molar density liquid", "molar_density_liquid");
  dependencies_.insert(KeyTag{dens_key_, tag});

  frac_cond_key_ = Keys::readKey(plist_, domain, "fractional conductance", "fractional_conductance");
  dependencies_.insert(KeyTag{frac_cond_key_, tag});

  drag_exp_key_ = Keys::readKey(plist_, domain, "drag exponent", "drag_exponent");
  dependencies_.insert(KeyTag{drag_exp_key_, tag});

  // create the model
  Teuchos::ParameterList& sublist = plist_.sublist("overland conductivity model");
  model_ = Teuchos::rcp(new ManningConductivityModel(sublist));
}


Teuchos::RCP<Evaluator>
OverlandConductivitySubgridEvaluator::Clone() const {
  return Teuchos::rcp(new OverlandConductivitySubgridEvaluator(*this));
}


// Required methods from EvaluatorSecondaryMonotypeCV
void OverlandConductivitySubgridEvaluator::Evaluate_(const State& S,
        const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> mobile_depth = S.GetPtr<CompositeVector>(mobile_depth_key_, tag);
  Teuchos::RCP<const CompositeVector> slope = S.GetPtr<CompositeVector>(slope_key_, tag);
  Teuchos::RCP<const CompositeVector> coef = S.GetPtr<CompositeVector>(coef_key_, tag);

  Teuchos::RCP<const CompositeVector> frac_cond = S.GetPtr<CompositeVector>(frac_cond_key_, tag);
  Teuchos::RCP<const CompositeVector> drag = S.GetPtr<CompositeVector>(drag_exp_key_, tag);

  Teuchos::RCP<const CompositeVector> dens = S.GetPtr<CompositeVector>(dens_key_, tag);

  std::string comp = "cell";
  const Epetra_MultiVector& mobile_depth_v = *mobile_depth->ViewComponent(comp,false);
  const Epetra_MultiVector& slope_v = *slope->ViewComponent(comp,false);
  const Epetra_MultiVector& coef_v = *coef->ViewComponent(comp,false);

  const Epetra_MultiVector& frac_cond_v = *frac_cond->ViewComponent(comp,false);
  const Epetra_MultiVector& drag_v = *drag->ViewComponent(comp,false);
  const Epetra_MultiVector& dens_v = *dens->ViewComponent(comp,false);
  Epetra_MultiVector& result_v = *result[0]->ViewComponent(comp,false);

  int ncomp = result[0]->size(comp, false);
  for (int i=0; i!=ncomp; ++i) {
    result_v[0][i] = model_->Conductivity(mobile_depth_v[0][i], slope_v[0][i], coef_v[0][i]);
    result_v[0][i] *= dens_v[0][i] * std::pow(frac_cond_v[0][i], drag_v[0][i] + 1);
  }

  // NOTE: Should deal with boundary face here! --etc
}

void OverlandConductivitySubgridEvaluator::EvaluatePartialDerivative_(
    const State& S,
    const Key& wrt_key, const Tag& wrt_tag, const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> mobile_depth = S.GetPtr<CompositeVector>(mobile_depth_key_, tag);
  Teuchos::RCP<const CompositeVector> slope = S.GetPtr<CompositeVector>(slope_key_, tag);
  Teuchos::RCP<const CompositeVector> coef = S.GetPtr<CompositeVector>(coef_key_, tag);

  Teuchos::RCP<const CompositeVector> frac_cond = S.GetPtr<CompositeVector>(frac_cond_key_, tag);
  Teuchos::RCP<const CompositeVector> drag = S.GetPtr<CompositeVector>(drag_exp_key_, tag);

  Teuchos::RCP<const CompositeVector> dens = S.GetPtr<CompositeVector>(dens_key_, tag);

  if (wrt_key == mobile_depth_key_) {
    std::string comp = "cell";
    const Epetra_MultiVector& mobile_depth_v = *mobile_depth->ViewComponent(comp,false);
    const Epetra_MultiVector& slope_v = *slope->ViewComponent(comp,false);
    const Epetra_MultiVector& coef_v = *coef->ViewComponent(comp,false);

    const Epetra_MultiVector& frac_cond_v = *frac_cond->ViewComponent(comp,false);
    const Epetra_MultiVector& drag_v = *drag->ViewComponent(comp,false);
    const Epetra_MultiVector& dens_v = *dens->ViewComponent(comp,false);
    Epetra_MultiVector& result_v = *result[0]->ViewComponent(comp,false);

    int ncomp = result[0]->size(comp, false);
    for (int i=0; i!=ncomp; ++i) {
      result_v[0][i] = model_->DConductivityDDepth(mobile_depth_v[0][i], slope_v[0][i], coef_v[0][i]);
      result_v[0][i] *= dens_v[0][i] * std::pow(frac_cond_v[0][i], drag_v[0][i] + 1);
    }

  } else if (wrt_key == dens_key_) {
    std::string comp = "cell";
    const Epetra_MultiVector& mobile_depth_v = *mobile_depth->ViewComponent(comp,false);
    const Epetra_MultiVector& slope_v = *slope->ViewComponent(comp,false);
    const Epetra_MultiVector& coef_v = *coef->ViewComponent(comp,false);

    const Epetra_MultiVector& frac_cond_v = *frac_cond->ViewComponent(comp,false);
    const Epetra_MultiVector& drag_v = *drag->ViewComponent(comp,false);
    const Epetra_MultiVector& dens_v = *dens->ViewComponent(comp,false);
    Epetra_MultiVector& result_v = *result[0]->ViewComponent(comp,false);

    int ncomp = result[0]->size(comp, false);
    for (int i=0; i!=ncomp; ++i) {
      result_v[0][i] = model_->Conductivity(mobile_depth_v[0][i], slope_v[0][i], coef_v[0][i]);
      result_v[0][i] *= std::pow(frac_cond_v[0][i], drag_v[0][i] + 1);
    }

  } else if (wrt_key == frac_cond_key_) {
    std::string comp = "cell";
    const Epetra_MultiVector& mobile_depth_v = *mobile_depth->ViewComponent(comp,false);
    const Epetra_MultiVector& slope_v = *slope->ViewComponent(comp,false);
    const Epetra_MultiVector& coef_v = *coef->ViewComponent(comp,false);

    const Epetra_MultiVector& frac_cond_v = *frac_cond->ViewComponent(comp,false);
    const Epetra_MultiVector& drag_v = *drag->ViewComponent(comp,false);
    const Epetra_MultiVector& dens_v = *dens->ViewComponent(comp,false);
    Epetra_MultiVector& result_v = *result[0]->ViewComponent(comp,false);

    int ncomp = result[0]->size(comp, false);
    for (int i=0; i!=ncomp; ++i) {
      result_v[0][i] = model_->Conductivity(mobile_depth_v[0][i], slope_v[0][i], coef_v[0][i]);
      result_v[0][i] *= dens_v[0][i] * (drag_v[0][i] + 1) *
        std::pow(frac_cond_v[0][i], drag_v[0][i]);
    }
  } else {
    result[0]->PutScalar(0.);
  }
}


void OverlandConductivitySubgridEvaluator::EnsureCompatibility_ToDeps_(State& S)
{
  auto akeytag = my_keys_[0];
  const auto& my_fac = S.Require<CompositeVector,CompositeVectorSpace>(akeytag.first, akeytag.second);
  if (my_fac.Mesh() != Teuchos::null) {
    CompositeVectorSpace dep_fac;
    dep_fac.SetMesh(my_fac.Mesh())
      ->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, 1);
    EvaluatorSecondaryMonotypeCV::EnsureCompatibility_ToDeps_(S, dep_fac);
  }
}

} // namespace Flow
} // namespace Amanzi

