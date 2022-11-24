/*
  The ideal gas equation of state model is an algebraic model with dependencies.

  Generated via evaluator_generator with:

    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_GENERAL_EOS_IDEAL_GAS_MODEL_HH_
#define AMANZI_GENERAL_EOS_IDEAL_GAS_MODEL_HH_

namespace Amanzi {
namespace General {
namespace Relations {

class EosIdealGasModel {
 public:
  explicit EosIdealGasModel(Teuchos::ParameterList& plist);

  double Density(double temp, double pres) const;

  double DDensityDTemperature(double temp, double pres) const;
  double DDensityDPressure(double temp, double pres) const;

 protected:
  void InitializeFromPlist_(Teuchos::ParameterList& plist);

 protected:
  double cv_;
  double T0_;
};

} // namespace Relations
} // namespace General
} // namespace Amanzi

#endif
