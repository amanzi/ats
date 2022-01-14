/*
  The evaporative flux relaxation evaluator is an algebraic evaluator of a given model.

  Generated via evaluator_generator with:
    myKeyFirst = evaporative
    evalNameString = evaporative flux relaxation
    evalNameCaps = EVAPORATIVE_FLUX_RELAXATION
    namespaceCaps = SURFACEBALANCE
    evalClassName = EvaporativeFluxRelaxation
    namespace = SurfaceBalance
    myMethodDeclarationArgs = double wc, double rho, double L
    myKey = evaporative_flux
    evalName = evaporative_flux_relaxation
    modelMethodDeclaration =   double EvaporativeFlux(double wc, double rho, double L) const;
    myKeyMethod = EvaporativeFlux
    myMethodArgs = wc_v[0][i], rho_v[0][i], L_v[0][i]
    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_SURFACEBALANCE_EVAPORATIVE_FLUX_RELAXATION_EVALUATOR_HH_
#define AMANZI_SURFACEBALANCE_EVAPORATIVE_FLUX_RELAXATION_EVALUATOR_HH_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class EvaporativeFluxRelaxationModel;

class EvaporativeFluxRelaxationEvaluator : public EvaluatorSecondaryMonotypeCV {

 public:
  explicit
  EvaporativeFluxRelaxationEvaluator(Teuchos::ParameterList& plist);
  EvaporativeFluxRelaxationEvaluator(const EvaporativeFluxRelaxationEvaluator& other);

  virtual Teuchos::RCP<Evaluator> Clone() const;

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  Teuchos::RCP<EvaporativeFluxRelaxationModel> get_model() { return model_; }

 protected:
  void InitializeFromPlist_();

  Key wc_key_;
  Key rho_key_;
  Key L_key_;

  Teuchos::RCP<EvaporativeFluxRelaxationModel> model_;

 private:
  static Utils::RegisteredFactory<Evaluator,EvaporativeFluxRelaxationEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace

#endif
