/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  The evaporative flux relaxation model is an algebraic model with dependencies.

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

*/

#ifndef AMANZI_SURFACEBALANCE_EVAPORATIVE_FLUX_RELAXATION_MODEL_HH_
#define AMANZI_SURFACEBALANCE_EVAPORATIVE_FLUX_RELAXATION_MODEL_HH_

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class EvaporativeFluxRelaxationModel {
 public:
  explicit EvaporativeFluxRelaxationModel(Teuchos::ParameterList& plist);

  double EvaporativeFlux(double wc, double rho, double L) const;

  double DEvaporativeFluxDLitterWaterContent(double wc, double rho, double L) const;
  double DEvaporativeFluxDSurfaceMolarDensityLiquid(double wc, double rho, double L) const;
  double DEvaporativeFluxDLitterThickness(double wc, double rho, double L) const;

 protected:
  void InitializeFromPlist_(Teuchos::ParameterList& plist);

 protected:
  double wc_sat_;
  double tau_;
};

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi

#endif
