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

#include "Teuchos_ParameterList.hpp"
#include "dbc.hh"
#include "latent_heat_model.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

// Constructor from ParameterList
LatentHeatModel::LatentHeatModel(Teuchos::ParameterList& plist)
{
  InitializeFromPlist_(plist);
}


// Initialize parameters
void
LatentHeatModel::InitializeFromPlist_(Teuchos::ParameterList& plist)
{
  Le_ = plist.get<double>("latent heat of vaporization [MJ/mol]", 0.0449994810744);
}


// main method
double
LatentHeatModel::LatentHeat(double qe) const
{
  AMANZI_ASSERT(0);
  return 0.;
}

double
LatentHeatModel::DLatentHeatDEvaporativeFlux(double qe) const
{
  AMANZI_ASSERT(0);
  return 0.;
}

} //namespace
} //namespace
} //namespace
  
