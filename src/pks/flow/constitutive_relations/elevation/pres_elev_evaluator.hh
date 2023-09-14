/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Evaluates the potential surface upon which overland flow acts.
/*!

.. math::
   h + z

`"evaluator type`" = 

.. _pres-elev-evaluator-spec:
.. admonition:: pres-elev-evaluator-spec

   KEYS:

   - `"height`" **DOMAIN-ponded_depth** Names the height variable. [m]
   - `"elevation`" **DOMAIN-elevation** Names the elevation variable. [m]


NOTE: This is a legacy evaluator, and is not in the factory, so need not be in
the input spec.  However, we include it here because this could easily be
abstracted for new potential surfaces, kinematic wave, etc, at which point it
would need to be added to the factory and the input spec.

NOTE: This could easily be replaced by a generic Additive_ Evaluator.

*/

#ifndef AMANZI_FLOWRELATIONS_PRES_ELEV_EVALUATOR_
#define AMANZI_FLOWRELATIONS_PRES_ELEV_EVALUATOR_

#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Flow {

class PresElevEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit PresElevEvaluator(Teuchos::ParameterList& plist);
  PresElevEvaluator(const PresElevEvaluator& other) = default;
  Teuchos::RCP<Evaluator> Clone() const override;

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

 private:
  Key pres_key_;
  Key elev_key_;
};

} // namespace Flow
} // namespace Amanzi

#endif
