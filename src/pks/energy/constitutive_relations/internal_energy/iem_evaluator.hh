/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  The IEM Evaluator simply calls the IEM with the correct arguments.

*/

/*!

Computes (specific) internal energy of as a function of temperature.

`"evaluator type`" = `"iem`"

.. _iem_evaluator-spec:
.. admonition:: iem_evaluator-spec

   * `"IEM parameters`" ``[IEM_model-typedinline-spec-list]``

   KEYS:

   - `"temperature`"

*/

#ifndef AMANZI_ENERGY_RELATIONS_IEM_EVALUATOR_
#define AMANZI_ENERGY_RELATIONS_IEM_EVALUATOR_

#include "Factory.hh"
#include "iem.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Energy {

class IEMEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  // constructor format for all derived classes
  explicit IEMEvaluator(Teuchos::ParameterList& plist);
  IEMEvaluator(Teuchos::ParameterList& plist, const Teuchos::RCP<IEM>& iem);
  IEMEvaluator(const IEMEvaluator& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override;

  Teuchos::RCP<IEM> get_IEM() { return iem_; }

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

  void InitializeFromPlist_();

  Key temp_key_;
  Teuchos::RCP<IEM> iem_;

 private:
  static Utils::RegisteredFactory<Evaluator, IEMEvaluator> factory_;
};

} // namespace Energy
} // namespace Amanzi

#endif
