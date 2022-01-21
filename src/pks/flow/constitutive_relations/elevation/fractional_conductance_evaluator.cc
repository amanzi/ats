/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

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
  const Epetra_MultiVector& vpd = *S.Get<CompositeVector>(vpd_key_, tag).ViewComponent("cell",false);
  const Epetra_MultiVector& del_max = *S.Get<CompositeVector>(delta_max_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& del_ex = *S.Get<CompositeVector>(delta_ex_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& depr_depth = *S.Get<CompositeVector>(depr_depth_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& mobile_depth = *S.Get<CompositeVector>(mobile_depth_key_, tag).ViewComponent("cell",false);
  Epetra_MultiVector& res = *result[0]->ViewComponent("cell",false);

  int ncells = res.MyLength();
  for (int c=0; c!=ncells; ++c) {
    if (mobile_depth[0][c] <= 0.0) {
      res[0][c] = 0;
    } else {
      double vol_depr_depth = Microtopography::volumetricDepth(depr_depth[0][c], del_max[0][c], del_ex[0][c]);
      res[0][c] = (vpd[0][c] - vol_depr_depth) / (mobile_depth[0][c]);
    }
  }
}


void
FractionalConductanceEvaluator::EvaluatePartialDerivative_(const State& S,
        const Key& wrt_key, const Tag& wrt_tag, const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  const Epetra_MultiVector& vpd = *S.Get<CompositeVector>(vpd_key_, tag).ViewComponent("cell",false);
  const Epetra_MultiVector& del_max = *S.Get<CompositeVector>(delta_max_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& del_ex = *S.Get<CompositeVector>(delta_ex_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& depr_depth = *S.Get<CompositeVector>(depr_depth_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& mobile_depth = *S.Get<CompositeVector>(mobile_depth_key_, tag).ViewComponent("cell",false);
  Epetra_MultiVector& res = *result[0]->ViewComponent("cell",false);

  if (wrt_key == mobile_depth_key_) {
    int ncells = res.MyLength();
    for (int c=0; c!=ncells; ++c) {
      if (mobile_depth[0][c] <= 0.0) {
        res[0][c] = 0;
      } else {
        double vol_depr_depth = Microtopography::volumetricDepth(depr_depth[0][c], del_max[0][c], del_ex[0][c]);
        res[0][c] = - (vpd[0][c] - vol_depr_depth) * std::pow(mobile_depth[0][c],-2.);
      }
    }

  } else if (wrt_key == vpd_key_) {
    int ncells = res.MyLength();
    for (int c=0; c!=ncells; ++c) {
      if (mobile_depth[0][c] <= 0.0) {
        res[0][c] = 0;
      } else {
        res[0][c] =  1.0 / mobile_depth[0][c];
      }
    }

  } else {
    Errors::Message msg("VolumetricHeightSubgridEvaluator: Not Implemented: no derivatives implemented other than ponded depth.");
    Exceptions::amanzi_throw(msg);
  }
}

} //namespace
} //namespace
} //namespace
