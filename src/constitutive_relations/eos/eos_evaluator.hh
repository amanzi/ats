/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  EOSFieldEvaluator is the interface between state/data and the model, an EOS.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_RELATIONS_EOS_EVALUATOR_HH_
#define AMANZI_RELATIONS_EOS_EVALUATOR_HH_

#include "eos.hh"
#include "factory.hh"
//#include "secondary_variables_field_evaluator.hh"
#include "EvaluatorSecondaries.hh"

namespace Amanzi {
namespace Relations {

class EOSEvaluator : public EvaluatorSecondaries {

 public:
  enum EOSMode { EOS_MODE_MASS, EOS_MODE_MOLAR, EOS_MODE_BOTH };

  // constructor format for all derived classes
  explicit
  EOSEvaluator(Teuchos::ParameterList& plist);

  EOSEvaluator(const EOSEvaluator& other);
  virtual Teuchos::RCP<Evaluator> Clone()  const override;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void Update_(State& S) override;
  virtual void UpdateDerivative_(State& S, const Key& wrt_key, const Key& wrt_tag) override;

  virtual bool IsDifferentiableWRT(const State& S, const Key& wrt_key, const Key& wrt_tag) const override {
    return false;
  }
  
  virtual void EnsureCompatibility(State& S) override {};
  
  Teuchos::RCP<EOS> get_EOS() { return eos_; }
 protected:
  // the actual model
  Teuchos::RCP<EOS> eos_;
  EOSMode mode_;

  Key tag_;
  // Keys for fields
  // dependencies
  Key temp_key_;
  Key pres_key_;

 private:
  static Utils::RegisteredFactory<Evaluator,EOSEvaluator> factory_;
};

} // namespace
} // namespace

#endif
