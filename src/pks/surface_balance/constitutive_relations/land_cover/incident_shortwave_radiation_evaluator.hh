/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
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
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class IncidentShortwaveRadiationModel;

class IncidentShortwaveRadiationEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit IncidentShortwaveRadiationEvaluator(Teuchos::ParameterList& plist);
  IncidentShortwaveRadiationEvaluator(const IncidentShortwaveRadiationEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;


  Teuchos::RCP<IncidentShortwaveRadiationModel> get_model() { return model_; }

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;
  void InitializeFromPlist_();

 protected:
  Key slope_key_;
  Key aspect_key_;
  Key qSWin_key_;

  Teuchos::RCP<IncidentShortwaveRadiationModel> model_;

 private:
  static Utils::RegisteredFactory<Evaluator, IncidentShortwaveRadiationEvaluator> reg_;
};

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
