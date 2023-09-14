/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Drainage rate from the canopy to the lower layers.
/*!

A simple model based on relaxation from current water content to a saturated water content.

.. code::
   
          |
          | source
          V
         /   \
      I /     \
       V       |
   --Theta--    | T
       ^       |
       | D     |
       V       V
   -- -- -- -- -- -- --


This is the model for drainage D.

Drainage is given by:

.. math::
   D = max(0, \frac{(\Theta - \Theta_sat)}{\tau})

.. _drainage-evaluator-spec:
.. admonition:: drainage-evaluator-spec

   * `"drainage timescale [s]`" ``[double]`` **864** Timescale over which drainage occurs.
   * `"saturated specific water content [m^3 H2O / m^2 leaf area]`" ``[double]`` **1e-4**
      The thickness of the wetting surface -- determines maximum canopy water storage.\

   KEYS:
   - "area index"
   - "water equivalent"

*/

#pragma once

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class DrainageEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  // constructor format for all derived classes
  explicit DrainageEvaluator(Teuchos::ParameterList& plist);

  DrainageEvaluator(const DrainageEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

 protected:
  Key drainage_key_;
  Key fracwet_key_;

  Key ai_key_;
  Key wc_key_;

  double tau_;
  double wc_sat_;
  double n_liq_;

 private:
  static Amanzi::Utils::RegisteredFactory<Evaluator, DrainageEvaluator> reg_;
};

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
