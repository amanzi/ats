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

#include "Teuchos_ParameterList.hpp"
#include "dbc.hh"
#include "evaporative_flux_relaxation_model.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

// Constructor from ParameterList
EvaporativeFluxRelaxationModel::EvaporativeFluxRelaxationModel(Teuchos::ParameterList& plist)
{
  InitializeFromPlist_(plist);
}


// Initialize parameters
void
EvaporativeFluxRelaxationModel::InitializeFromPlist_(Teuchos::ParameterList& plist)
{
  wc_sat_ = plist.get<double>("saturated litter water content [-]");
  tau_ = plist.get<double>("drying time [s]");
}


// main method
double
EvaporativeFluxRelaxationModel::EvaporativeFlux(double wc, double rho, double L) const
{
  return wc / (wc_sat_ * L * rho) / tau_;
}

double
EvaporativeFluxRelaxationModel::DEvaporativeFluxDLitterWaterContent(double wc,
                                                                    double rho,
                                                                    double L) const
{
  return 1.0 / (wc_sat_ * L * rho) / tau_;
}

double
EvaporativeFluxRelaxationModel::DEvaporativeFluxDSurfaceMolarDensityLiquid(double wc,
                                                                           double rho,
                                                                           double L) const
{
  return -wc / (wc_sat_ * L * rho) / tau_ / rho;
}

double
EvaporativeFluxRelaxationModel::DEvaporativeFluxDLitterThickness(double wc,
                                                                 double rho,
                                                                 double L) const
{
  return -wc / (wc_sat_ * L * rho) / tau_ / L;
}

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
