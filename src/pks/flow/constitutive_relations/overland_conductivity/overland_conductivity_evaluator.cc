/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluates the conductivity of surface flow.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "overland_conductivity_evaluator.hh"
#include "manning_conductivity_model.hh"
#include "split_denominator_conductivity_model.hh"
#include "ponded_depth_passthrough_conductivity_model.hh"

namespace Amanzi {
namespace Flow {

OverlandConductivityEvaluator::OverlandConductivityEvaluator(Teuchos::ParameterList& plist)
    : EvaluatorSecondaryMonotypeCV(plist)
{
  Key domain = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  if (plist_.isParameter("height key") || plist_.isParameter("ponded depth key")
      || plist_.isParameter("height key suffix") || plist_.isParameter("ponded depth key suffix")) {
    Errors::Message message("OverlandConductivity: only use \"depth key\" or \"depth key suffix\", not \"height key\" or \"ponded depth key\".");
    Exceptions::amanzi_throw(message);
  }
  depth_key_ = Keys::readKey(plist_, domain, "depth", "ponded_depth");
  dependencies_.insert(KeyTag{depth_key_, tag});

  slope_key_ = Keys::readKey(plist_, domain, "slope", "slope_magnitude");
  dependencies_.insert(KeyTag{slope_key_, tag});

  coef_key_ = Keys::readKey(plist_, domain, "coefficient", "manning_coefficient");
  dependencies_.insert(KeyTag{coef_key_, tag});

  dt_swe_factor_ = plist_.get<double>("dt factor [s]", -1);
  if (dt_swe_factor_ > 0) {
    double swe_factor = plist_.get<double>("swe density factor [-]", 10.0);
    dt_swe_factor_ *= swe_factor;
  }

  dens_ = plist_.get<bool>("include density", true);
  if (dens_) {
    dens_key_ = Keys::readKey(plist_, domain, "molar density liquid", "molar_density_liquid");
    dependencies_.insert(KeyTag{dens_key_, tag});
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
void OverlandConductivityEvaluator::Evaluate_(const State& S,
        const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> depth = S.GetPtr<CompositeVector>(depth_key_, tag);
  Teuchos::RCP<const CompositeVector> slope = S.GetPtr<CompositeVector>(slope_key_, tag);
  Teuchos::RCP<const CompositeVector> coef = S.GetPtr<CompositeVector>(coef_key_, tag);

#ifdef ENABLE_DBC
  double min_coef = 1.;
  coef->MinValue(&min_coef);
  if (min_coef <= 1.e-12) {
    Errors::Message message("Overland Conductivity Evaluator: Manning coeficient has at least one value that is non-positive.  Perhaps you forgot to set the \"boundary_face\" component?");
    Exceptions::amanzi_throw(message);
  }
#endif

  for (const auto& comp : *result[0]) {
    const Epetra_MultiVector& depth_v = *depth->ViewComponent(comp,false);
    const Epetra_MultiVector& slope_v = *slope->ViewComponent(comp,false);
    const Epetra_MultiVector& coef_v = *coef->ViewComponent(comp,false);
    Epetra_MultiVector& result_v = *result[0]->ViewComponent(comp,false);

    int ncomp = result[0]->size(comp, false);
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
      const Epetra_MultiVector& dens_v = *S.Get<CompositeVector>(dens_key_, tag).ViewComponent(comp,false);
      for (int i=0; i!=ncomp; ++i) result_v[0][i] *= dens_v[0][i];
    }
  }
}


void
OverlandConductivityEvaluator::EvaluatePartialDerivative_(const State& S,
        const Key& wrt_key, const Tag& wrt_tag, const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> depth = S.GetPtr<CompositeVector>(depth_key_, tag);
  Teuchos::RCP<const CompositeVector> slope = S.GetPtr<CompositeVector>(slope_key_, tag);
  Teuchos::RCP<const CompositeVector> coef = S.GetPtr<CompositeVector>(coef_key_, tag);

  if (wrt_key == depth_key_) {
    for (const auto& comp : *result[0]) {
      const Epetra_MultiVector& depth_v = *depth->ViewComponent(comp,false);
      const Epetra_MultiVector& slope_v = *slope->ViewComponent(comp,false);
      const Epetra_MultiVector& coef_v = *coef->ViewComponent(comp,false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(comp,false);

      int ncomp = result[0]->size(comp, false);
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
        const Epetra_MultiVector& dens_v = *S.Get<CompositeVector>(dens_key_, tag).ViewComponent(comp,false);
        for (int i=0; i!=ncomp; ++i) {
          result_v[0][i] *= dens_v[0][i];
        }
      }
    }

  } else if (wrt_key == dens_key_) {
    AMANZI_ASSERT(dens_);
    for (const auto& comp : *result[0]) {
      const Epetra_MultiVector& depth_v = *depth->ViewComponent(comp,false);
      const Epetra_MultiVector& slope_v = *slope->ViewComponent(comp,false);
      const Epetra_MultiVector& coef_v = *coef->ViewComponent(comp,false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(comp,false);

      int ncomp = result[0]->size(comp, false);
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

  } else {
    // FIX ME -- need to add derivatives of conductivity model wrt slope, coef --etc
    result[0]->PutScalar(0.);
  }
}


} //namespace
} //namespace

