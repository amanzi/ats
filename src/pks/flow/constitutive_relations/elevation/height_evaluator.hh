/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
/*!

Computes ponded depth from surface water pressure.

.. math::
   h = \frac{H(p - p_{atm})}{\rho g}

where :math:`H` is the Heaviside function.

`"evaluator type`" = `"ponded depth`"

.. _evaluator-ponded-depth-spec:
.. admonition:: evaluator-ponded-depth-spec

   KEYS:

   - `"pressure`"
   - `"mass density`"
   - `"atmospheric pressure`"
   - `"gravity`"
*/

#ifndef AMANZI_FLOW_RELATIONS_HEIGHT_EVALUATOR_
#define AMANZI_FLOW_RELATIONS_HEIGHT_EVALUATOR_

#include "EvaluatorSecondaryMonotype.hh"
#include "Factory.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {

class HeightModel;

class HeightEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  // constructor format for all derived classes
  explicit HeightEvaluator(Teuchos::ParameterList& plist);
  HeightEvaluator(const HeightEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  Teuchos::RCP<HeightModel> get_Model() { return model_; }

  void set_bar(bool bar) { bar_ = bar; }

 protected:
  // Needs a special EnsureCompatibility to get around trying to find face
  // values and derivatives of face values.
  virtual void EnsureCompatibility_ToDeps_(State& S) override;

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

 protected:
  Key dens_key_;
  Key pres_key_;
  Key gravity_key_;
  Key patm_key_;
  bool bar_;

  Teuchos::RCP<HeightModel> model_;

 private:
  static Utils::RegisteredFactory<Evaluator, HeightEvaluator> factory_;
};

} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi

#endif
