/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  EOSFieldEvaluator is the interface between state/data and the model, an EOS.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_RELATIONS_EOS_EVALUATOR_HH_
#define AMANZI_RELATIONS_EOS_EVALUATOR_HH_

#include "eos.hh"
#include "Factory.hh"
#include "secondary_variables_field_evaluator.hh"

namespace Amanzi {
namespace Relations {

class EOSEvaluator : public SecondaryVariablesFieldEvaluator {

 public:
  enum EOSMode { EOS_MODE_MASS, EOS_MODE_MOLAR, EOS_MODE_BOTH };

  // constructor format for all derived classes
  explicit
  EOSEvaluator(Teuchos::ParameterList& plist);

  EOSEvaluator(const EOSEvaluator& other);
  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const std::vector<Teuchos::Ptr<CompositeVector> >& results);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
                                               Key wrt_key, const std::vector<Teuchos::Ptr<CompositeVector> >& results);

  Teuchos::RCP<EOS> get_EOS() { return eos_; }
 protected:
  // the actual model
  Teuchos::RCP<EOS> eos_;
  EOSMode mode_;
};

} // namespace
} // namespace

#endif
