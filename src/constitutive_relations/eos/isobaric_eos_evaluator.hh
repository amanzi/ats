/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  EOSEvaluator is the interface between state/data and the model, an EOS.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_RELATIONS_ISOBARIC_EOS_EVALUATOR_HH_
#define AMANZI_RELATIONS_ISOBARIC_EOS_EVALUATOR_HH_

#include "eos.hh"
#include "factory.hh"
//#include "secondary_variables_field_evaluator.hh"
#include "EvaluatorSecondaries.hh"

namespace Amanzi {
namespace Relations {

class IsobaricEOSEvaluator : public EvaluatorSecondaries {

 public:
  enum EOSMode { EOS_MODE_MASS, EOS_MODE_MOLAR, EOS_MODE_BOTH };

  // constructor format for all derived classes
  explicit
  IsobaricEOSEvaluator(Teuchos::ParameterList& plist);

  IsobaricEOSEvaluator(const IsobaricEOSEvaluator& other);
  virtual Teuchos::RCP<Evaluator> Clone() const;

  Teuchos::RCP<EOS> get_EOS() { return eos_; }

protected:
  // the actual model
  Teuchos::RCP<EOS> eos_;
  EOSMode mode_;

  // Required methods from SecondaryVariableEvaluator
  // These do the actual work
  virtual void Evaluate_(const State& S,
                         const std::vector<Teuchos::Ptr<CompositeVector> > & results);
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key, const Key& wrt_tag,
                                          const std::vector<Teuchos::Ptr<CompositeVector> > & results);
  
  // Keys for fields
  // dependencies
  Key temp_key_;
  Key pres_key_;
  Key a_key_;

 private:
  static Utils::RegisteredFactory<Evaluator,IsobaricEOSEvaluator> factory_;
};

} // namespace
} // namespace

#endif
