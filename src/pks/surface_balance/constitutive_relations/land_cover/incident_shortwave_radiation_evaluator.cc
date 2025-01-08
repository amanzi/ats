/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

//! Evaluates shortwave as a function of slope/aspect/etc.
/*!

Aspect modified shortwave radiation is determined by a factor which
is multiplied by the 'incoming radiation incident on a flat surface'
to determine the 'incoming radiation incident on a sloping surface of
a given aspect' as a function of latitude, slope, aspect, and Julian
day of the year, and time of day.

Note that some careful checking and experimentation has found that, in
general, the daily average incident radiation times the 12-noon aspect
modifier correlates reasonably well with the daily average of the
product of the hourly incident radiation and the hourly aspect
modifier.  It is notably better than the daily average radiation times
the daily average aspect modifier.

Derived from LandLab code, which is released under the MIT license:
https://github.com/landlab/landlab/blob/master/landlab/components/radiation/radiation.py

.. _incident_shortwave_radiation_evaluator-spec:
.. admonition:: incident_shortwave_radiation_evaluator-spec

    * `"incident shortwave radiation parameters`" ``[incident_shortwave_radiation_model-spec]``

    KEYS:
    * `"slope`"
    * `"aspect`"
    * `"incoming shortwave radiation`"

*/

#include "incident_shortwave_radiation_evaluator.hh"
#include "incident_shortwave_radiation_model.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

// Constructor from ParameterList
IncidentShortwaveRadiationEvaluator::IncidentShortwaveRadiationEvaluator(
  Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  Teuchos::ParameterList& sublist = plist_.sublist("incident shortwave radiation parameters");
  model_ = Teuchos::rcp(new IncidentShortwaveRadiationModel(sublist));
  InitializeFromPlist_();
}


// Virtual copy constructor
Teuchos::RCP<Evaluator>
IncidentShortwaveRadiationEvaluator::Clone() const
{
  return Teuchos::rcp(new IncidentShortwaveRadiationEvaluator(*this));
}


// Initialize by setting up dependencies
void
IncidentShortwaveRadiationEvaluator::InitializeFromPlist_()
{
  // Set up my dependencies
  // - defaults to prefixed via domain
  Tag tag = my_keys_.front().second;
  Key domain_name = Keys::getDomain(my_keys_.front().first);

  // - pull Keys from plist
  // dependency: slope
  slope_key_ = Keys::readKey(plist_, domain_name, "slope magnitude", "slope_magnitude");
  dependencies_.insert(KeyTag{ slope_key_, tag });

  // dependency: aspect
  aspect_key_ = Keys::readKey(plist_, domain_name, "aspect", "aspect");
  dependencies_.insert(KeyTag{ aspect_key_, tag });

  // dependency: incoming_shortwave_radiation
  qSWin_key_ = Keys::readKey(
    plist_, domain_name, "incoming shortwave radiation", "incoming_shortwave_radiation");
  dependencies_.insert(KeyTag{ qSWin_key_, tag });
}


void
IncidentShortwaveRadiationEvaluator::Evaluate_(const State& S,
                                               const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> slope = S.GetPtr<CompositeVector>(slope_key_, tag);
  Teuchos::RCP<const CompositeVector> aspect = S.GetPtr<CompositeVector>(aspect_key_, tag);
  Teuchos::RCP<const CompositeVector> qSWin = S.GetPtr<CompositeVector>(qSWin_key_, tag);

  for (CompositeVector::name_iterator comp = result[0]->begin(); comp != result[0]->end(); ++comp) {
    const Epetra_MultiVector& slope_v = *slope->ViewComponent(*comp, false);
    const Epetra_MultiVector& aspect_v = *aspect->ViewComponent(*comp, false);
    const Epetra_MultiVector& qSWin_v = *qSWin->ViewComponent(*comp, false);
    Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);
    double time = S.get_time();

    int ncomp = result[0]->size(*comp, false);
    for (int i = 0; i != ncomp; ++i) {
      result_v[0][i] =
        model_->IncidentShortwaveRadiation(slope_v[0][i], aspect_v[0][i], qSWin_v[0][i], time);
    }
  }
}


void
IncidentShortwaveRadiationEvaluator::EvaluatePartialDerivative_(
  const State& S,
  const Key& wrt_key,
  const Tag& wrt_tag,
  const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> slope = S.GetPtr<CompositeVector>(slope_key_, tag);
  Teuchos::RCP<const CompositeVector> aspect = S.GetPtr<CompositeVector>(aspect_key_, tag);
  Teuchos::RCP<const CompositeVector> qSWin = S.GetPtr<CompositeVector>(qSWin_key_, tag);
  double time = S.get_time();

  if (wrt_key == slope_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin(); comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& slope_v = *slope->ViewComponent(*comp, false);
      const Epetra_MultiVector& aspect_v = *aspect->ViewComponent(*comp, false);
      const Epetra_MultiVector& qSWin_v = *qSWin->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = model_->DIncidentShortwaveRadiationDSlope(
          slope_v[0][i], aspect_v[0][i], qSWin_v[0][i], time);
      }
    }

  } else if (wrt_key == aspect_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin(); comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& slope_v = *slope->ViewComponent(*comp, false);
      const Epetra_MultiVector& aspect_v = *aspect->ViewComponent(*comp, false);
      const Epetra_MultiVector& qSWin_v = *qSWin->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = model_->DIncidentShortwaveRadiationDAspect(
          slope_v[0][i], aspect_v[0][i], qSWin_v[0][i], time);
      }
    }

  } else if (wrt_key == qSWin_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin(); comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& slope_v = *slope->ViewComponent(*comp, false);
      const Epetra_MultiVector& aspect_v = *aspect->ViewComponent(*comp, false);
      const Epetra_MultiVector& qSWin_v = *qSWin->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = model_->DIncidentShortwaveRadiationDIncomingShortwaveRadiation(
          slope_v[0][i], aspect_v[0][i], qSWin_v[0][i], time);
      }
    }

  } else {
    AMANZI_ASSERT(false);
  }
}


} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
