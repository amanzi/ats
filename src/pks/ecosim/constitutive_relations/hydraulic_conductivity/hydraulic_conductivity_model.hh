/*
  The hydraulic conductivity model is an algebraic model with dependencies.

  Generated via evaluator_generator with:
Richards water content evaluator: the standard form as a function of liquid saturation.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_ECOSIM_HYDRAULIC_CONDUCTIVITY_MODEL_HH_
#define AMANZI_ECOSIM_HYDRAULIC_CONDUCTIVITY_MODEL_HH_

namespace Amanzi {
namespace Ecosim {
namespace Relations {

class HydraulicConductivityModel {

 public:
  explicit
  HydraulicConductivityModel(Teuchos::ParameterList& plist);

  double HydraulicConductivity(double k, double rho, double mu, double gz) const;

  double DHydraulicConductivityDPermeability(double k, double rho, double mu, double gz) const;
  double DHydraulicConductivityDMassDensityLiquid(double k, double rho, double mu, double gz) const;
  double DHydraulicConductivityDViscosityLiquid(double k, double rho, double mu, double gz) const;

 protected:
  void InitializeFromPlist_(Teuchos::ParameterList& plist);

 protected:

  double g_;

};

} //namespace
} //namespace
} //namespace

#endif
