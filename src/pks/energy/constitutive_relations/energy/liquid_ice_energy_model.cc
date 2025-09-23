/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  The liquid+ice energy model is an algebraic model with dependencies.

  Generated via evaluator_generator with:
Energy for a two-phase, liquid+ice evaluator.

*/

#include "Teuchos_ParameterList.hpp"
#include "dbc.hh"
#include "liquid_ice_energy_model.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Energy {
namespace Relations {

// Constructor from ParameterList
LiquidIceEnergyModel::LiquidIceEnergyModel(Teuchos::ParameterList& plist)
{
  InitializeFromPlist_(plist);
}


// Initialize parameters
void
LiquidIceEnergyModel::InitializeFromPlist_(Teuchos::ParameterList& plist)
{}


// main method
double
LiquidIceEnergyModel::Energy(double phi,
                             double phi0,
                             double sl,
                             double nl,
                             double ul,
                             double si,
                             double ni,
                             double ui,
                             double rho_r,
                             double ur,
                             double cv) const
{
  return cv * (phi * (ni * si * ui + nl * sl * ul) + rho_r * ur * (-phi0 + 1));
}

double
LiquidIceEnergyModel::DEnergyDPorosity(double phi,
                                       double phi0,
                                       double sl,
                                       double nl,
                                       double ul,
                                       double si,
                                       double ni,
                                       double ui,
                                       double rho_r,
                                       double ur,
                                       double cv) const
{
  return cv * (ni * si * ui + nl * sl * ul);
}

double
LiquidIceEnergyModel::DEnergyDBasePorosity(double phi,
                                           double phi0,
                                           double sl,
                                           double nl,
                                           double ul,
                                           double si,
                                           double ni,
                                           double ui,
                                           double rho_r,
                                           double ur,
                                           double cv) const
{
  return -cv * rho_r * ur;
}

double
LiquidIceEnergyModel::DEnergyDSaturationLiquid(double phi,
                                               double phi0,
                                               double sl,
                                               double nl,
                                               double ul,
                                               double si,
                                               double ni,
                                               double ui,
                                               double rho_r,
                                               double ur,
                                               double cv) const
{
  return cv * nl * phi * ul;
}

double
LiquidIceEnergyModel::DEnergyDMolarDensityLiquid(double phi,
                                                 double phi0,
                                                 double sl,
                                                 double nl,
                                                 double ul,
                                                 double si,
                                                 double ni,
                                                 double ui,
                                                 double rho_r,
                                                 double ur,
                                                 double cv) const
{
  return cv * phi * sl * ul;
}

double
LiquidIceEnergyModel::DEnergyDInternalEnergyLiquid(double phi,
                                                   double phi0,
                                                   double sl,
                                                   double nl,
                                                   double ul,
                                                   double si,
                                                   double ni,
                                                   double ui,
                                                   double rho_r,
                                                   double ur,
                                                   double cv) const
{
  return cv * nl * phi * sl;
}

double
LiquidIceEnergyModel::DEnergyDSaturationIce(double phi,
                                            double phi0,
                                            double sl,
                                            double nl,
                                            double ul,
                                            double si,
                                            double ni,
                                            double ui,
                                            double rho_r,
                                            double ur,
                                            double cv) const
{
  return cv * ni * phi * ui;
}

double
LiquidIceEnergyModel::DEnergyDMolarDensityIce(double phi,
                                              double phi0,
                                              double sl,
                                              double nl,
                                              double ul,
                                              double si,
                                              double ni,
                                              double ui,
                                              double rho_r,
                                              double ur,
                                              double cv) const
{
  return cv * phi * si * ui;
}

double
LiquidIceEnergyModel::DEnergyDInternalEnergyIce(double phi,
                                                double phi0,
                                                double sl,
                                                double nl,
                                                double ul,
                                                double si,
                                                double ni,
                                                double ui,
                                                double rho_r,
                                                double ur,
                                                double cv) const
{
  return cv * ni * phi * si;
}

double
LiquidIceEnergyModel::DEnergyDDensityRock(double phi,
                                          double phi0,
                                          double sl,
                                          double nl,
                                          double ul,
                                          double si,
                                          double ni,
                                          double ui,
                                          double rho_r,
                                          double ur,
                                          double cv) const
{
  return cv * ur * (-phi0 + 1);
}

double
LiquidIceEnergyModel::DEnergyDInternalEnergyRock(double phi,
                                                 double phi0,
                                                 double sl,
                                                 double nl,
                                                 double ul,
                                                 double si,
                                                 double ni,
                                                 double ui,
                                                 double rho_r,
                                                 double ur,
                                                 double cv) const
{
  return cv * rho_r * (-phi0 + 1);
}

double
LiquidIceEnergyModel::DEnergyDCellVolume(double phi,
                                         double phi0,
                                         double sl,
                                         double nl,
                                         double ul,
                                         double si,
                                         double ni,
                                         double ui,
                                         double rho_r,
                                         double ur,
                                         double cv) const
{
  return phi * (ni * si * ui + nl * sl * ul) + rho_r * ur * (-phi0 + 1);
}

} // namespace Relations
} // namespace Energy
} // namespace ATS_Physics
} // namespace Amanzi
