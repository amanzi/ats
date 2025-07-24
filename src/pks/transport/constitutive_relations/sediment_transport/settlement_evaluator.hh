/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatsky (dasvyat@lanl.gov)
*/

/*!

The settlement evaluator gets the settlement rates. We estimate the deposition
due to settling, Q_ds which is mainly due to the formation and breakup of
flocs, by way of the formulation proposed by Einstein and Krone [1962] who
assumed that most of the particles settled in flocs as long as near bed shear
stresses are small enough to prevent their breaking up, namely,

.. math::

   Q_ds = \left\{\begin{array}{ll}
                    \rho_s w_s min(C, 0.5)\left(1-\frac{\tau_0}{\tau_d}\right) &\qquad\mbox{if} \qquad \tau_0<\tau_d \\
                    0  &\qquad\mbox{if} \qquad \tau_0\ge\tau_d
                 \end{array} \right.

where :math:`\rho_s` is the sediment density, :math:`w_s` is the settling velocity which depends on the size
of sediment flocs, :math:`\tau_d` is a critical shear stress, C in a sediment molar ratio.

`"evaluator type`" = `"settlement rate`"

.. _evaluator-settlement-rate-spec:
.. admonition:: evaluator-settlement-rate-spec

   * `"critical shear stress`" ``[double]`` **0.1**
   * `"settling velocity`" ``[double]`` **0.0001**
   * `"drag coefficient`" ``[double]`` **0.02**
   * `"specific weight of water`" ``[double]`` **9806.0**

   DEPENDENCIES:

   - `"velocity`" **DOMAIN-velocity**
   - `"sediment`" **DOMAIN-sediment**

*/

#ifndef AMANZI_SETTLEMENTRATE_EVALUATOR_
#define AMANZI_SETTLEMENTRATE_EVALUATOR_

// TPLs
#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {

class SettlementRateEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit SettlementRateEvaluator(Teuchos::ParameterList& plist);
  SettlementRateEvaluator(const SettlementRateEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

  double tau_d_;
  double ws_;
  double gamma_;
  double lambda_, umax_, xi_;
  double sediment_density_;
  double Cf_;
  Key velocity_key_, sediment_key_;

  static Utils::RegisteredFactory<Evaluator, SettlementRateEvaluator> factory_;
};

} // namespace Amanzi

#endif
