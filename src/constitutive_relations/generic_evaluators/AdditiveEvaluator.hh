/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! A generic evaluator for summing a collection of fields.
/*!

.. _evaluator-additive-spec:
.. admonition:: evaluator-additive-spec
   * `"constant shift`" ``[double]`` **0** A constant value to add to the sum.

   * `"enforce positivity`" ``[bool]`` **false** If true, max the result with 0.

   * `"DEPENDENCY coefficient`" ``[double]`` **1.0** A multiple for each dependency in
     the list below.

   ONE OF

   * `"dependencies`" ``[Array(string)]`` The things to sum.

   OR

   * `"evaluator dependency suffixes`" ``[Array(string)]``

   END

*/

#pragma once

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Relations {

class AdditiveEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  // constructor format for all derived classes
  explicit AdditiveEvaluator(Teuchos::ParameterList& plist);

  AdditiveEvaluator(const AdditiveEvaluator& other) = default;
  Teuchos::RCP<Evaluator> Clone() const override;

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  void EvaluatePartialDerivative_(const State& S,
                                  const Key& wrt_key,
                                  const Tag& wrt_tag,
                                  const std::vector<CompositeVector*>& result) override;

 protected:
  std::map<Key, double> coefs_;
  double shift_;
  bool positive_;

 private:
  static Utils::RegisteredFactory<Evaluator, AdditiveEvaluator> factory_;
};

} // namespace Relations
} // namespace Amanzi
