/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ahmad Jan (jana@ornl.gov)
*/

//! Evaluates snow melt via USDA - Natural Resources Conservation Service model
/*!

From:  National Engineering Handbook (NEH) part 630, Chapter 11

Snow melt rate is given by:

.. math::
   SM = H(T_{snow}^{expected} - 273.15) R

where :math:`R` is the snow melt rate per degree-day and
:math:`T_{snow}^{expected}` is the expected snow temperature, which is
typically given by :math:`T_{snow}^{expected} = T_{air} - \Delta`, where
:math:`\Delta` is the expected air-snow temperature difference [C] (note this
is NOT prescribed here -- the user must supply the expected snow temperature
via an evalutor).

Note that the Heaviside function is used to ensure this is only active when the
expected snow temperature is above 0 C.

Then, a linear transition factor is applied to ensure that this goes to zero as
the snow SWE goes to zero.  That factor is 1 at the snow transition depth, and
0 when snow SWE is 0.  This uses LandCover for the snow_ground_transition
parameter.

.. _snow-meltrate-evaluator-spec:
.. admonition:: snow-meltrate-evaluator-spec

   * `"snow melt rate [mm day^-1 C^-1]`" ``[double]`` **2.74**
     the melt rate per degree-day, above 0 C, e.g. :math:`R` above.

   KEYS:

   - `"snow water equivalent`" **DOMAIN-water_equivalent**
   - `"expected snow temperature`"  **DOMAIN-expected_temperature**

*/

#pragma once

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "LandCover.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class SnowMeltRateEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit SnowMeltRateEvaluator(Teuchos::ParameterList& plist);
  SnowMeltRateEvaluator(const SnowMeltRateEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new SnowMeltRateEvaluator(*this));
  }

 protected:
  virtual void EnsureCompatibility_ToDeps_(State& S) override;

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

 protected:
  Key exp_temp_key_;
  Key snow_key_;

  double melt_rate_;
  double snow_temp_shift_;

  Key domain_;
  bool compatibility_checked_;

  LandCoverMap land_cover_;

 private:
  static Utils::RegisteredFactory<Evaluator, SnowMeltRateEvaluator> reg_;
};

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
