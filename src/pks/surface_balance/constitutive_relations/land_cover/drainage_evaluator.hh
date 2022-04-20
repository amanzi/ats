/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! Drainage rate from the canopy to the lower layers.

/*!

A simple model based on relaxation from current water content to a saturated water content.

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
  ----------------------

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
#include "secondary_variables_field_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class DrainageEvaluator : public SecondaryVariablesFieldEvaluator {

 public:

  // constructor format for all derived classes
  explicit
  DrainageEvaluator(Teuchos::ParameterList& plist);

  DrainageEvaluator(const DrainageEvaluator& other) = default;
  virtual Teuchos::RCP<FieldEvaluator> Clone() const override;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const std::vector<Teuchos::Ptr<CompositeVector> >& results) override;
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const std::vector<Teuchos::Ptr<CompositeVector> > & results) override;

 protected:
  Key drainage_key_;
  Key fracwet_key_;

  Key ai_key_;
  Key wc_key_;

  double tau_;
  double wc_sat_;
  double n_liq_;

 private:
  static Amanzi::Utils::RegisteredFactory<FieldEvaluator,DrainageEvaluator> reg_;

};

} // namespace
} // namespace
} // namespace

