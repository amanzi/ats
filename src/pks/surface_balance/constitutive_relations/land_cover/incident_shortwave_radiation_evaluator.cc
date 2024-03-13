/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
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

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

const std::string IncidentShortwaveRadiationEvaluator::eval_type = "incident shortwave radiation";

// Constructor from ParameterList
IncidentShortwaveRadiationEvaluator::IncidentShortwaveRadiationEvaluator(
  const Teuchos::RCP<Teuchos::ParameterList>& plist)
  : EvaluatorModelCV<IncidentShortwaveRadiationModel>(plist)
{
  doy0_ = plist->sublist("model parameters").get<int>("day of year at time 0 [Julian days]", 0);
  if (doy0_ < 0 || doy0_ > 364) {
    Errors::Message msg("IncidentShortwaveRadiationModel: \"day of year at time 0 [Julian days]\" "
                        "not in valid range [0,364]");
    Exceptions::amanzi_throw(msg);
  }
}


// Virtual copy constructor
Teuchos::RCP<Evaluator>
IncidentShortwaveRadiationEvaluator::Clone() const
{
  return Teuchos::rcp(new IncidentShortwaveRadiationEvaluator(*this));
}


void
IncidentShortwaveRadiationEvaluator::Evaluate_(const State& S,
                                               const std::vector<CompositeVector*>& result)
{
  double time_days = S.get_time(my_keys_.front().second) / 86400;
  double doy = std::fmod((double)doy0_ + time_days, (double)365);
  int doy_i = std::lround(doy);
  if (doy_i == 365) {
    // can round up!
    doy_i = 0;
    doy = doy - 365.0;
  }

  model_->doy = doy;
  model_->doy_i = doy_i;
  EvaluatorModelCV<IncidentShortwaveRadiationModel>::Evaluate_(S, result);
}


} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
