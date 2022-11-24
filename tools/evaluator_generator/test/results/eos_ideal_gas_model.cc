/*
  The ideal gas equation of state model is an algebraic model with dependencies.

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
    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "Teuchos_ParameterList.hpp"
#include "dbc.hh"
#include "eos_ideal_gas_model.hh"

namespace Amanzi {
namespace General {
namespace Relations {

// Constructor from ParameterList
EosIdealGasModel::EosIdealGasModel(Teuchos::ParameterList& plist)
{
  InitializeFromPlist_(plist);
}


// Initialize parameters
void
EosIdealGasModel::InitializeFromPlist_(Teuchos::ParameterList& plist)
{
  cv_ = plist.get<double>("heat capacity");
}


// main method
double
EosIdealGasModel::Density(double temp, double pres) const
{
  return AMANZI_ASSERT(False);
}

double
EosIdealGasModel::DDensityDTemperature(double temp, double pres) const
{
  return AMANZI_ASSERT(False);
}

double
EosIdealGasModel::DDensityDPressure(double temp, double pres) const
{
  return AMANZI_ASSERT(False);
}

} // namespace Relations
} // namespace General
} // namespace Amanzi
