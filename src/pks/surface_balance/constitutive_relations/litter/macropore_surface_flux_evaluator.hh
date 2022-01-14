/*
  The macropore-surface flux evaluator is an algebraic evaluator of a given model.

  Generated via evaluator_generator with:
    evalName = macropore_surface_flux
    modelMethodDeclaration =   double MacroporeSurfaceFlux(double pM, double ps, double krs, double krM, double K) const;
    namespaceCaps = SURFACEBALANCE
    namespace = SurfaceBalance
    evalNameCaps = MACROPORE_SURFACE_FLUX
    myMethodArgs = pM_v[0][i], ps_v[0][i], krs_v[0][i], krM_v[0][i], K_v[0][i]
    myKeyMethod = MacroporeSurfaceFlux
    myKeyFirst = macropore
    evalNameString = macropore-surface flux
    myMethodDeclarationArgs = double pM, double ps, double krs, double krM, double K
    evalClassName = MacroporeSurfaceFlux
    myKey = macropore_surface_flux
    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_SURFACEBALANCE_MACROPORE_SURFACE_FLUX_EVALUATOR_HH_
#define AMANZI_SURFACEBALANCE_MACROPORE_SURFACE_FLUX_EVALUATOR_HH_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class MacroporeSurfaceFluxModel;

class MacroporeSurfaceFluxEvaluator : public EvaluatorSecondaryMonotypeCV {

 public:
  explicit
  MacroporeSurfaceFluxEvaluator(Teuchos::ParameterList& plist);
  MacroporeSurfaceFluxEvaluator(const MacroporeSurfaceFluxEvaluator& other);

  virtual Teuchos::RCP<Evaluator> Clone() const;

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  Teuchos::RCP<MacroporeSurfaceFluxModel> get_model() { return model_; }

 protected:
  void InitializeFromPlist_();

  Key pM_key_;
  Key ps_key_;
  Key krs_key_;
  Key krM_key_;
  Key K_key_;

  Teuchos::RCP<MacroporeSurfaceFluxModel> model_;

 private:
  static Utils::RegisteredFactory<Evaluator,MacroporeSurfaceFluxEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace

#endif
