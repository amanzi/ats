/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*!

A non-isothermal viscosity model intended for use within a range of
temperatures from well below freezing to ~100C.

.. _viscosity-evaluator-spec:
.. admonition:: viscosity-evaluator-spec

   * `"viscosity model parameters`" ``[viscosity-typedinline-spec-list]``

   KEYS:

   - `"temperature`"

 */


#ifndef AMANZI_RELATIONS_VISC_EVALUATOR_HH_
#define AMANZI_RELATIONS_VISC_EVALUATOR_HH_

#include "viscosity_relation.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Relations {

class ViscosityEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  // constructor format for all derived classes
  explicit ViscosityEvaluator(Teuchos::ParameterList& plist);

  ViscosityEvaluator(const ViscosityEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

 protected:
  // the actual model
  Teuchos::RCP<ViscosityRelation> visc_;

  // Keys for fields
  // dependencies
  Key temp_key_;

 private:
  static Utils::RegisteredFactory<Evaluator, ViscosityEvaluator> factory_;
};

} // namespace Relations
} // namespace Amanzi

#endif
