/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/
//! Evaluates shortwave as a function of slope/aspect/etc.

/*!
  Generated via evaluator_generator with:
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
*/

#pragma once

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class IncidentShortwaveRadiationModel;

class IncidentShortwaveRadiationEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  IncidentShortwaveRadiationEvaluator(Teuchos::ParameterList& plist);
  IncidentShortwaveRadiationEvaluator(const IncidentShortwaveRadiationEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  Teuchos::RCP<IncidentShortwaveRadiationModel> get_model() { return model_; }

 protected:
  void InitializeFromPlist_();

  Key slope_key_;
  Key aspect_key_;
  Key qSWin_key_;

  Teuchos::RCP<IncidentShortwaveRadiationModel> model_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,IncidentShortwaveRadiationEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace

