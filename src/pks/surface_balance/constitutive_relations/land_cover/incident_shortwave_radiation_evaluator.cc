/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
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
IncidentShortwaveRadiationEvaluator::IncidentShortwaveRadiationEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist)
{
  Teuchos::ParameterList& sublist = plist_.sublist("incident shortwave radiation parameters");
  model_ = Teuchos::rcp(new IncidentShortwaveRadiationModel(sublist));
  InitializeFromPlist_();
}


// Copy constructor
IncidentShortwaveRadiationEvaluator::IncidentShortwaveRadiationEvaluator(const IncidentShortwaveRadiationEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    slope_key_(other.slope_key_),
    aspect_key_(other.aspect_key_),
    qSWin_key_(other.qSWin_key_),    
    model_(other.model_) {}


// Virtual copy constructor
Teuchos::RCP<FieldEvaluator>
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
  Key domain_name = Keys::getDomain(my_key_);

  // - pull Keys from plist
  // dependency: slope
  slope_key_ = Keys::readKey(plist_, domain_name, "slope magnitude", "slope_magnitude");
  dependencies_.insert(slope_key_);

  // dependency: aspect
  aspect_key_ = Keys::readKey(plist_, domain_name, "aspect", "aspect");
  dependencies_.insert(aspect_key_);

  // dependency: incoming_shortwave_radiation
  qSWin_key_ = Keys::readKey(plist_, domain_name, "incoming shortwave radiation", "incoming_shortwave_radiation");
  dependencies_.insert(qSWin_key_);
}


void
IncidentShortwaveRadiationEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
  Teuchos::RCP<const CompositeVector> slope = S->GetFieldData(slope_key_);
  Teuchos::RCP<const CompositeVector> aspect = S->GetFieldData(aspect_key_);
  Teuchos::RCP<const CompositeVector> qSWin = S->GetFieldData(qSWin_key_);

  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    const Epetra_MultiVector& slope_v = *slope->ViewComponent(*comp, false);
    const Epetra_MultiVector& aspect_v = *aspect->ViewComponent(*comp, false);
    const Epetra_MultiVector& qSWin_v = *qSWin->ViewComponent(*comp, false);
    Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);
    double time = S->time();

    int ncomp = result->size(*comp, false);
    for (int i=0; i!=ncomp; ++i) {
      result_v[0][i] = model_->IncidentShortwaveRadiation(slope_v[0][i], aspect_v[0][i], qSWin_v[0][i], time);
    }
  }
}


void
IncidentShortwaveRadiationEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  Teuchos::RCP<const CompositeVector> slope = S->GetFieldData(slope_key_);
  Teuchos::RCP<const CompositeVector> aspect = S->GetFieldData(aspect_key_);
  Teuchos::RCP<const CompositeVector> qSWin = S->GetFieldData(qSWin_key_);
  double time = S->time();

  if (wrt_key == slope_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& slope_v = *slope->ViewComponent(*comp, false);
      const Epetra_MultiVector& aspect_v = *aspect->ViewComponent(*comp, false);
      const Epetra_MultiVector& qSWin_v = *qSWin->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DIncidentShortwaveRadiationDSlope(slope_v[0][i], aspect_v[0][i], qSWin_v[0][i], time);
      }
    }

  } else if (wrt_key == aspect_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& slope_v = *slope->ViewComponent(*comp, false);
      const Epetra_MultiVector& aspect_v = *aspect->ViewComponent(*comp, false);
      const Epetra_MultiVector& qSWin_v = *qSWin->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DIncidentShortwaveRadiationDAspect(slope_v[0][i], aspect_v[0][i], qSWin_v[0][i], time);
      }
    }

  } else if (wrt_key == qSWin_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& slope_v = *slope->ViewComponent(*comp, false);
      const Epetra_MultiVector& aspect_v = *aspect->ViewComponent(*comp, false);
      const Epetra_MultiVector& qSWin_v = *qSWin->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DIncidentShortwaveRadiationDIncomingShortwaveRadiation(slope_v[0][i], aspect_v[0][i], qSWin_v[0][i], time);
      }
    }

  } else {
    AMANZI_ASSERT(false);
  }
}


} //namespace
} //namespace
} //namespace
