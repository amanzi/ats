/*
  The latent heat from evaporative flux model is an algebraic model with dependencies.

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

#ifndef AMANZI_SURFACEBALANCE_LATENT_HEAT_MODEL_HH_
#define AMANZI_SURFACEBALANCE_LATENT_HEAT_MODEL_HH_

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class LatentHeatModel {

 public:
  explicit
  LatentHeatModel(Teuchos::ParameterList& plist);

  double LatentHeat(double qe) const;

  double DLatentHeatDEvaporativeFlux(double qe) const;
  
 protected:
  void InitializeFromPlist_(Teuchos::ParameterList& plist);

 protected:

  double Le_;

};

} //namespace
} //namespace
} //namespace

#endif
