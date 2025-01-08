/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

//! Fraction of incoming water that is intercepted.
#include "interception_fraction_evaluator.hh"
#include "interception_fraction_model.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

// Constructor from ParameterList
InterceptionFractionEvaluator::InterceptionFractionEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  Teuchos::ParameterList& sublist = plist_.sublist("interception fraction parameters");
  model_ = Teuchos::rcp(new InterceptionFractionModel(sublist));
  InitializeFromPlist_();
}


// Virtual copy constructor
Teuchos::RCP<Evaluator>
InterceptionFractionEvaluator::Clone() const
{
  return Teuchos::rcp(new InterceptionFractionEvaluator(*this));
}


// Initialize by setting up dependencies
void
InterceptionFractionEvaluator::InitializeFromPlist_()
{
  Key akey = my_keys_.front().first;
  Tag tag = my_keys_.front().second;
  Key domain = Keys::getDomain(akey);
  akey = Keys::getVarName(akey);
  Key domain_surf = Keys::readDomainHint(plist_, domain, "canopy", "surface");
  Key domain_snow = Keys::readDomainHint(plist_, domain, "canopy", "snow");
  my_keys_.clear();

  // my keys
  interception_key_ = Keys::in(akey, "interception") ? akey : "interception";
  interception_key_ = Keys::readKey(plist_, domain, "interception", interception_key_);
  my_keys_.emplace_back(KeyTag{ interception_key_, tag });

  throughfall_rain_key_ =
    (Keys::in(akey, "rain") && Keys::in(akey, "throughfall")) ? akey : "throughfall_drainage_rain";
  throughfall_rain_key_ =
    Keys::readKey(plist_, domain, "throughfall and drainage rain", throughfall_rain_key_);
  my_keys_.emplace_back(KeyTag{ throughfall_rain_key_, tag });

  throughfall_snow_key_ =
    (Keys::in(akey, "snow") && Keys::in(akey, "throughfall")) ? akey : "throughfall_drainage_snow";
  throughfall_snow_key_ =
    Keys::readKey(plist_, domain, "throughfall and drainage snow", throughfall_snow_key_);
  my_keys_.emplace_back(KeyTag{ throughfall_snow_key_, tag });

  // - pull Keys from plist
  // dependency: surface-area_index
  ai_key_ = Keys::readKey(plist_, domain, "area index", "area_index");
  dependencies_.insert(KeyTag{ ai_key_, tag });
  rain_key_ = Keys::readKey(plist_, domain_surf, "precipitation rain", "precipitation_rain");
  dependencies_.insert(KeyTag{ rain_key_, tag });
  snow_key_ = Keys::readKey(plist_, domain_snow, "precipitation snow", "precipitation");
  dependencies_.insert(KeyTag{ snow_key_, tag });
  drainage_key_ = Keys::readKey(plist_, domain, "drainage", "drainage");
  dependencies_.insert(KeyTag{ drainage_key_, tag });
  air_temp_key_ = Keys::readKey(plist_, domain_surf, "air temperature", "air_temperature");
  dependencies_.insert(KeyTag{ air_temp_key_, tag });
}


void
InterceptionFractionEvaluator::Evaluate_(const State& S,
                                         const std::vector<CompositeVector*>& results)
{
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> ai = S.GetPtr<CompositeVector>(ai_key_, tag);
  Teuchos::RCP<const CompositeVector> rain = S.GetPtr<CompositeVector>(rain_key_, tag);
  Teuchos::RCP<const CompositeVector> snow = S.GetPtr<CompositeVector>(snow_key_, tag);
  Teuchos::RCP<const CompositeVector> drainage = S.GetPtr<CompositeVector>(drainage_key_, tag);
  Teuchos::RCP<const CompositeVector> air_temp = S.GetPtr<CompositeVector>(air_temp_key_, tag);

  for (CompositeVector::name_iterator comp = results[0]->begin(); comp != results[0]->end();
       ++comp) {
    const Epetra_MultiVector& ai_v = *ai->ViewComponent(*comp, false);
    const Epetra_MultiVector& rain_v = *rain->ViewComponent(*comp, false);
    const Epetra_MultiVector& snow_v = *snow->ViewComponent(*comp, false);
    const Epetra_MultiVector& drainage_v = *drainage->ViewComponent(*comp, false);
    const Epetra_MultiVector& air_temp_v = *air_temp->ViewComponent(*comp, false);

    Epetra_MultiVector& inter_v = *results[0]->ViewComponent(*comp, false);
    Epetra_MultiVector& tfr_v = *results[1]->ViewComponent(*comp, false);
    Epetra_MultiVector& tfs_v = *results[2]->ViewComponent(*comp, false);

    for (int i = 0; i != inter_v.MyLength(); ++i) {
      double coef = model_->InterceptionFraction(ai_v[0][i]);
      double total_precip = rain_v[0][i] + snow_v[0][i];
      inter_v[0][i] = total_precip * coef;

      double frac_r =
        total_precip > 0 ? rain_v[0][i] / total_precip : (air_temp_v[0][i] > 273.15 ? 1 : 0);
      tfr_v[0][i] = (1 - coef) * rain_v[0][i] + frac_r * drainage_v[0][i];
      tfs_v[0][i] = (1 - coef) * snow_v[0][i] + (1 - frac_r) * drainage_v[0][i];
    }
  }
}


void
InterceptionFractionEvaluator::EvaluatePartialDerivative_(
  const State& S,
  const Key& wrt_key,
  const Tag& wrt_tag,
  const std::vector<CompositeVector*>& results)
{
  if (wrt_key == drainage_key_) {
    Tag tag = my_keys_.front().second;
    Teuchos::RCP<const CompositeVector> ai = S.GetPtr<CompositeVector>(ai_key_, tag);
    Teuchos::RCP<const CompositeVector> rain = S.GetPtr<CompositeVector>(rain_key_, tag);
    Teuchos::RCP<const CompositeVector> snow = S.GetPtr<CompositeVector>(snow_key_, tag);
    Teuchos::RCP<const CompositeVector> drainage = S.GetPtr<CompositeVector>(drainage_key_, tag);
    Teuchos::RCP<const CompositeVector> air_temp = S.GetPtr<CompositeVector>(air_temp_key_, tag);

    for (CompositeVector::name_iterator comp = results[0]->begin(); comp != results[0]->end();
         ++comp) {
      const Epetra_MultiVector& ai_v = *ai->ViewComponent(*comp, false);
      const Epetra_MultiVector& rain_v = *rain->ViewComponent(*comp, false);
      const Epetra_MultiVector& snow_v = *snow->ViewComponent(*comp, false);
      const Epetra_MultiVector& drainage_v = *drainage->ViewComponent(*comp, false);
      const Epetra_MultiVector& air_temp_v = *air_temp->ViewComponent(*comp, false);

      Epetra_MultiVector& inter_v = *results[0]->ViewComponent(*comp, false);
      Epetra_MultiVector& tfr_v = *results[1]->ViewComponent(*comp, false);
      Epetra_MultiVector& tfs_v = *results[2]->ViewComponent(*comp, false);

      for (int i = 0; i != inter_v.MyLength(); ++i) {
        inter_v[0][i] = 0.;

        double total_precip = rain_v[0][i] + snow_v[0][i];
        double frac_r =
          total_precip > 0 ? rain_v[0][i] / total_precip : (air_temp_v[0][i] > 273.15 ? 1 : 0);
        tfr_v[0][i] = frac_r;
        tfs_v[0][i] = 1 - frac_r;
      }
    }
  }
}


} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
