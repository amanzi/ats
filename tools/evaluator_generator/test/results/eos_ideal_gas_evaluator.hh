/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  The ideal gas equation of state evaluator is an algebraic evaluator of a given model.

  Generated via evaluator_generator with:
    modelInitializeParamsList =   cv_ = plist.get<double>("heat capacity");
    evalName = eos_ideal_gas
    modelMethodDeclaration =   double Density(double temp, double pres) const;
    namespaceCaps = GENERAL
    namespace = General
    paramDeclarationList =   double cv_;
    evalNameCaps = EOS_IDEAL_GAS
    myMethodArgs = temp_v[0][i], pres_v[0][i]
    myKeyMethod = Density
    myKeyFirst = density
    evalNameString = ideal gas equation of state
    myMethodDeclarationArgs = double temp, double pres
    evalClassName = EosIdealGas
    myKey = density

*/

#ifndef AMANZI_GENERAL_EOS_IDEAL_GAS_EVALUATOR_HH_
#define AMANZI_GENERAL_EOS_IDEAL_GAS_EVALUATOR_HH_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace General {
namespace Relations {

class EosIdealGasModel;

class EosIdealGasEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit EosIdealGasEvaluator(Teuchos::ParameterList& plist);
  EosIdealGasEvaluator(const EosIdealGasEvaluator& other);

  virtual Teuchos::RCP<Evaluator> Clone() const;

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
                              const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
                                               Key wrt_key,
                                               const Teuchos::Ptr<CompositeVector>& result);

  Teuchos::RCP<EosIdealGasModel> get_model() { return model_; }

 protected:
  void InitializeFromPlist_();

  Key temp_key_;
  Key pres_key_;

  Teuchos::RCP<EosIdealGasModel> model_;

 private:
  static Utils::RegisteredFactory<Evaluator, EosIdealGasEvaluator> reg_;
};

} // namespace Relations
} // namespace General
} // namespace Amanzi

#endif
