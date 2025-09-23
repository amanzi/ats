/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  The liquid+gas energy model is an algebraic model with dependencies.

  Generated via evaluator_generator with:
Energy for a two-phase, liquid+water vapor evaluator.

*/

#include "Teuchos_ParameterList.hpp"
#include "dbc.hh"
#include "liquid_gas_energy_model.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Energy {
namespace Relations {

// Constructor from ParameterList
LiquidGasEnergyModel::LiquidGasEnergyModel(Teuchos::ParameterList& plist)
{
  InitializeFromPlist_(plist);
}


// Initialize parameters
void
LiquidGasEnergyModel::InitializeFromPlist_(Teuchos::ParameterList& plist)
{}


// main method
double
LiquidGasEnergyModel::Energy(double phi,
                             double phi0,
                             double sl,
                             double nl,
                             double ul,
                             double sg,
                             double ng,
                             double ug,
                             double rho_r,
                             double ur,
                             double cv) const
{
  return cv * (phi * (ng * sg * ug + nl * sl * ul) + rho_r * ur * (-phi0 + 1));
}

double
LiquidGasEnergyModel::DEnergyDPorosity(double phi,
                                       double phi0,
                                       double sl,
                                       double nl,
                                       double ul,
                                       double sg,
                                       double ng,
                                       double ug,
                                       double rho_r,
                                       double ur,
                                       double cv) const
{
  return cv * (ng * sg * ug + nl * sl * ul);
}

double
LiquidGasEnergyModel::DEnergyDBasePorosity(double phi,
                                           double phi0,
                                           double sl,
                                           double nl,
                                           double ul,
                                           double sg,
                                           double ng,
                                           double ug,
                                           double rho_r,
                                           double ur,
                                           double cv) const
{
  return -cv * rho_r * ur;
}

double
LiquidGasEnergyModel::DEnergyDSaturationLiquid(double phi,
                                               double phi0,
                                               double sl,
                                               double nl,
                                               double ul,
                                               double sg,
                                               double ng,
                                               double ug,
                                               double rho_r,
                                               double ur,
                                               double cv) const
{
  return cv * nl * phi * ul;
}

double
LiquidGasEnergyModel::DEnergyDMolarDensityLiquid(double phi,
                                                 double phi0,
                                                 double sl,
                                                 double nl,
                                                 double ul,
                                                 double sg,
                                                 double ng,
                                                 double ug,
                                                 double rho_r,
                                                 double ur,
                                                 double cv) const
{
  return cv * phi * sl * ul;
}

double
LiquidGasEnergyModel::DEnergyDInternalEnergyLiquid(double phi,
                                                   double phi0,
                                                   double sl,
                                                   double nl,
                                                   double ul,
                                                   double sg,
                                                   double ng,
                                                   double ug,
                                                   double rho_r,
                                                   double ur,
                                                   double cv) const
{
  return cv * nl * phi * sl;
}

double
LiquidGasEnergyModel::DEnergyDSaturationGas(double phi,
                                            double phi0,
                                            double sl,
                                            double nl,
                                            double ul,
                                            double sg,
                                            double ng,
                                            double ug,
                                            double rho_r,
                                            double ur,
                                            double cv) const
{
  return cv * ng * phi * ug;
}

double
LiquidGasEnergyModel::DEnergyDMolarDensityGas(double phi,
                                              double phi0,
                                              double sl,
                                              double nl,
                                              double ul,
                                              double sg,
                                              double ng,
                                              double ug,
                                              double rho_r,
                                              double ur,
                                              double cv) const
{
  return cv * phi * sg * ug;
}

double
LiquidGasEnergyModel::DEnergyDInternalEnergyGas(double phi,
                                                double phi0,
                                                double sl,
                                                double nl,
                                                double ul,
                                                double sg,
                                                double ng,
                                                double ug,
                                                double rho_r,
                                                double ur,
                                                double cv) const
{
  return cv * ng * phi * sg;
}

double
LiquidGasEnergyModel::DEnergyDDensityRock(double phi,
                                          double phi0,
                                          double sl,
                                          double nl,
                                          double ul,
                                          double sg,
                                          double ng,
                                          double ug,
                                          double rho_r,
                                          double ur,
                                          double cv) const
{
  return cv * ur * (-phi0 + 1);
}

double
LiquidGasEnergyModel::DEnergyDInternalEnergyRock(double phi,
                                                 double phi0,
                                                 double sl,
                                                 double nl,
                                                 double ul,
                                                 double sg,
                                                 double ng,
                                                 double ug,
                                                 double rho_r,
                                                 double ur,
                                                 double cv) const
{
  return cv * rho_r * (-phi0 + 1);
}

double
LiquidGasEnergyModel::DEnergyDCellVolume(double phi,
                                         double phi0,
                                         double sl,
                                         double nl,
                                         double ul,
                                         double sg,
                                         double ng,
                                         double ug,
                                         double rho_r,
                                         double ur,
                                         double cv) const
{
  return phi * (ng * sg * ug + nl * sl * ul) + rho_r * ur * (-phi0 + 1);
}

} // namespace Relations
} // namespace Energy
} // namespace ATS_Physics
} // namespace Amanzi
