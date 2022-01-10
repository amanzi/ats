/*
  The latent heat from evaporative flux evaluator is an algebraic evaluator of a given model.

  Generated via evaluator_generator with:
    keyDeclarationList =   Key qe_key_;
    myKeyFirst = latent
    evalNameString = latent heat from evaporative flux
    keyCopyConstructorList =     qe_key_(other.qe_key_),
    evalNameCaps = LATENT_HEAT
    namespaceCaps = SURFACEBALANCE
    paramDeclarationList =   double Le_;
    modelDerivDeclarationList =   double DLatentHeatDEvaporativeFlux(double qe) const;
    evalClassName = LatentHeat
    keyCompositeVectorList =   Teuchos::RCP<const CompositeVector> qe = S->GetPtr<CompositeVector>("evaporative_flux");
    namespace = SurfaceBalance
    modelInitializeParamsList =   Le_ = plist.get<double>("latent heat of vaporization [MJ/mol]", 0.0449994810744);
    myMethodDeclarationArgs = double qe
    myKey = latent_heat
    evalName = latent_heat
    modelMethodDeclaration =   double LatentHeat(double qe) const;
    myKeyMethod = LatentHeat
    myMethodArgs = qe_v[0][i]
    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_SURFACEBALANCE_LATENT_HEAT_EVALUATOR_HH_
#define AMANZI_SURFACEBALANCE_LATENT_HEAT_EVALUATOR_HH_

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class LatentHeatModel;

class LatentHeatEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  LatentHeatEvaluator(Teuchos::ParameterList& plist);
  LatentHeatEvaluator(const LatentHeatEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  Teuchos::RCP<LatentHeatModel> get_model() { return model_; }

 protected:
  void InitializeFromPlist_();

  Key qe_key_;

  Teuchos::RCP<LatentHeatModel> model_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,LatentHeatEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace

#endif
