/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
/*!

Computes enthalpy `[MJ mol^-1]` of as a function of internal energy, pressure, and density.

.. math::

   e = u + 10^{-6} * \frac{p}{n_l}

`"evaluator type`" = `"enthalpy`"

.. _evaluator-enthalpy-spec:
.. admonition:: evaluator-enthalpy-spec

   * `"include work term`" ``[bool]`` **false** If false, e = u, ignoring the work term.

   KEYS:

   - `"internal energy`"
   - `"pressure`"
   - `"mass density`"

*/

#ifndef AMANZI_ENTHALPY_EVALUATOR_HH_
#define AMANZI_ENTHALPY_EVALUATOR_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Energy {

class EnthalpyEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit EnthalpyEvaluator(Teuchos::ParameterList& plist);
  EnthalpyEvaluator(const EnthalpyEvaluator& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override;

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

 protected:
  Key pres_key_;
  Key dens_key_;
  Key ie_key_;
  bool include_work_;

 private:
  static Utils::RegisteredFactory<Evaluator, EnthalpyEvaluator> factory_;
};

} // namespace Energy
} // namespace Amanzi

#endif
