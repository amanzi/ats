/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Evaluates the radiation incident on a non-flat surface.
/*!

Aspect modified shortwave radiation is determined by a factor which is
multiplied by the 'incoming radiation incident on a flat surface' to determine
the 'incoming radiation incident on a sloping surface of a given aspect' as a
function of slope and aspect, Julian day of the year, and time of day.  The
latitude and Julian day of the year are used to modify this with both time of
day and seasonal changes of the planet.

Note that some careful checking and experimentation has found that, in
general, the daily average incoming radiation times the 12-noon aspect
modifier correlates reasonably well with the daily average of the
product of the hourly incoming radiation and the hourly aspect
modifier.  It is notably better than the daily average radiation times
the daily average aspect modifier.

This implementation is derived from `LandLab code
<https://github.com/landlab/landlab/blob/master/landlab/components/radiation/radiation.py>`_,
which is released under the MIT license.


.. _incident_shortwave_radiation_evaluator-spec:
.. admonition:: incident_shortwave_radiation_evaluator-spec

    * `"incident shortwave radiation parameters`" ``[incident_shortwave_radiation_model-spec]``

    KEYS:
    - `"slope`" **DOMAIN-slope_magnitude**
    - `"aspect`" **DOMAIN-aspect**
    - `"incoming shortwave radiation`" **DOMAIN-incoming_shortwave_radiation**

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
