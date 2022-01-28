/*
  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! EOSEvaluator is the interface between state/data and the model, an EOS.

#ifndef AMANZI_RELATIONS_EOS_EVALUATOR_HH_
#define AMANZI_RELATIONS_EOS_EVALUATOR_HH_

#include "eos.hh"
#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Relations {

class EOSEvaluator : public EvaluatorSecondaryMonotypeCV {

 public:
  enum EOSMode { EOS_MODE_MASS, EOS_MODE_MOLAR, EOS_MODE_BOTH };

  // constructor format for all derived classes
  explicit EOSEvaluator(Teuchos::ParameterList& plist);
  EOSEvaluator(const EOSEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  Teuchos::RCP<EOS> get_EOS() { return eos_; }

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S,
          const std::vector<CompositeVector*>& results) override;
  virtual void EvaluatePartialDerivative_(const State& S,
          const Key& wrt_key, const Tag& wrt_tag,
          const std::vector<CompositeVector*>& results) override;

 protected:
  void ParsePlistKeys_();
  void ParsePlistTemp_();
  void ParsePlistPres_();
  void ParsePlistConc_();

 protected:
  // the actual model
  Teuchos::RCP<EOS> eos_;

  Key temp_key_, pres_key_, conc_key_;
  EOSMode mode_;
  bool updated_once_;

 private:
  static Utils::RegisteredFactory<Evaluator,EOSEvaluator> fac_;
};

} // namespace
} // namespace

#endif
