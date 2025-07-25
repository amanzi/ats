/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatsky (dasvyat@lanl.gov)
*/

/*!

The erosion evaluator gets the erosion rates.  We evaluate the erosion flux,
:math:`Q_e`, by way of a relationship which can be applied when the bed
properties are relatively uniform over the depth and the bed is consolidated
[e.g., Mehta, 1984], namely,

.. math::

   Q_e = \left\{\begin{array}{ll}
                   Q_{e_0}\left(\frac{\tau_0}{\tau_e} - 1\right) &\qquad\mbox{if} \qquad \tau_0>\tau_e \\
                   0  &\qquad\mbox{if} \qquad \tau_0\le\tau_e
                \end{array} \right.

where :math:`Q_{e_0}` is an empirical coefficient,  :math:`\tau_e` is a critical shear stress

`"evaluator type`" = `"erosion rate`"

.. _evaluator-erosion-rate-spec:
.. admonition:: evaluator-erosion-rate-spec

   * `"critical shear stress`" ``[double]`` **0.4**
   * `"empirical coefficient`" ``[double]`` **3.0e-4**
   * `"drag coefficient`" ``[double]`` **0.02**
   * `"specific weight of water`" ``[double]`` **9806.0**

   DEPENDENCIES:

   - `"velocity`" **DOMAIN-velocity**

*/

#ifndef AMANZI_EROSIONRATE_EVALUATOR_
#define AMANZI_EROSIONRATE_EVALUATOR_

// TPLs
#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {

class ErosionRateEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit ErosionRateEvaluator(Teuchos::ParameterList& plist);
  ErosionRateEvaluator(const ErosionRateEvaluator& other);
  Teuchos::RCP<Evaluator> Clone() const override;

  // virtual void EvaluateElevationAndSlope_(const Teuchos::Ptr<State>& S,
  //         const std::vector<Teuchos::Ptr<CompositeVector> >& results) = 0;
  // virtual bool HasFieldChanged(const Teuchos::Ptr<State>& S, Key request);
  // virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S){};

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

  double tau_e_;
  double Qe_0_;
  double gamma_;
  double lambda_, umax_, xi_;
  double Cf_;

  Key velocity_key_;

  static Utils::RegisteredFactory<Evaluator, ErosionRateEvaluator> factory_;
};

} // namespace Amanzi

#endif
