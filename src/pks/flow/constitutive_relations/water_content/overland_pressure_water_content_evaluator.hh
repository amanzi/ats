/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
/*!

Computes water content for a water-only surface water problem.  This can be
used to both compute an extensive water content (which is always non-negative)
or a "bar" water content which can be negative for pressures less than
atmospheric pressure for use in preconditioners.

By default:

.. math::

   \Theta_s = \begin{cases}
                0 & \text{if} \, p < p_{atm} \\
                V \frac{1}{M_{H2O} g_z} \frac{(p - p_{atm})^2}{2 R} & \text{if} \, p_{atm} <= p < R \\
                V \frac{1}{M_{H2O} g_z} (p - p_{atm} - \frac{R}{2}) & \text{if} \, \text{R <= p}
              \end{cases}

for molar mass :math:`M_{H2O}`, rollover cutoff :math:`R`, and gravity :math:`g`.

If `"allow negative water content`" is true, then the first conditional is
ignored (in this case :math:`R = 0` and is ignored).

Note this is not valid for concentration-dependent density problems.

`"evaluator type`" = `"overland pressure water content`"

.. _evaluator-overland-pressure-water-content-spec:
.. admonition:: evaluator-overland-pressure-water-content-spec

  * `"allow negative water content`" ``[bool]`` **false**  See above

  * `"water content rollover [Pa]`" ``[double]`` **0** A smoothing term that
    also can represent subgrid topography.  Not typically used.

  * `"molar mass [kg mol^-1]`" ``[double]`` **0.0180153**

  DEPENDENCIES:

  - `"pressure`"
  - `"cell volume`"
  - `"atmospheric_pressure`"
  - `"gravity`"

*/

#ifndef AMANZI_FLOW_RELATIONS_OVERLAND_HEAD_WATER_CONTENT_EVALUATOR_
#define AMANZI_FLOW_RELATIONS_OVERLAND_HEAD_WATER_CONTENT_EVALUATOR_

#include "EvaluatorSecondaryMonotype.hh"
#include "Factory.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {

class OverlandPressureWaterContentEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  // constructor format for all derived classes
  explicit OverlandPressureWaterContentEvaluator(Teuchos::ParameterList& plist);
  OverlandPressureWaterContentEvaluator(const OverlandPressureWaterContentEvaluator& other) =
    default;

  virtual Teuchos::RCP<Evaluator> Clone() const override;

 protected:
  void InitializeFromPlist_();

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

 protected:
  Key pres_key_, cv_key_;

  double M_;
  bool bar_; // bar'd variable indicates this is potentially negative for
             // pressures less than atmospheric
  double rollover_;

 private:
  static Utils::RegisteredFactory<Evaluator, OverlandPressureWaterContentEvaluator> reg_;
};

} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi

#endif
