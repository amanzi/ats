/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  The richards energy model is an algebraic model with dependencies.

  Generated via evaluator_generator with:
Richards energy: the standard form as a function of liquid saturation and specific internal energy.

*/

#include "Teuchos_ParameterList.hpp"
#include "dbc.hh"
#include "richards_energy_model.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Energy {
namespace Relations {

// Constructor from ParameterList
RichardsEnergyModel::RichardsEnergyModel(Teuchos::ParameterList& plist)
{
  InitializeFromPlist_(plist);
}


// Initialize parameters
void
RichardsEnergyModel::InitializeFromPlist_(Teuchos::ParameterList& plist)
{}


// main method
double
RichardsEnergyModel::Energy(double phi,
                            double phi0,
                            double sl,
                            double nl,
                            double ul,
                            double rho_r,
                            double ur,
                            double cv) const
{
  return cv * (nl * phi * sl * ul + rho_r * ur * (-phi0 + 1));
}

double
RichardsEnergyModel::DEnergyDPorosity(double phi,
                                      double phi0,
                                      double sl,
                                      double nl,
                                      double ul,
                                      double rho_r,
                                      double ur,
                                      double cv) const
{
  return cv * nl * sl * ul;
}

double
RichardsEnergyModel::DEnergyDBasePorosity(double phi,
                                          double phi0,
                                          double sl,
                                          double nl,
                                          double ul,
                                          double rho_r,
                                          double ur,
                                          double cv) const
{
  return -cv * rho_r * ur;
}

double
RichardsEnergyModel::DEnergyDSaturationLiquid(double phi,
                                              double phi0,
                                              double sl,
                                              double nl,
                                              double ul,
                                              double rho_r,
                                              double ur,
                                              double cv) const
{
  return cv * nl * phi * ul;
}

double
RichardsEnergyModel::DEnergyDMolarDensityLiquid(double phi,
                                                double phi0,
                                                double sl,
                                                double nl,
                                                double ul,
                                                double rho_r,
                                                double ur,
                                                double cv) const
{
  return cv * phi * sl * ul;
}

double
RichardsEnergyModel::DEnergyDInternalEnergyLiquid(double phi,
                                                  double phi0,
                                                  double sl,
                                                  double nl,
                                                  double ul,
                                                  double rho_r,
                                                  double ur,
                                                  double cv) const
{
  return cv * nl * phi * sl;
}

double
RichardsEnergyModel::DEnergyDDensityRock(double phi,
                                         double phi0,
                                         double sl,
                                         double nl,
                                         double ul,
                                         double rho_r,
                                         double ur,
                                         double cv) const
{
  return cv * ur * (-phi0 + 1);
}

double
RichardsEnergyModel::DEnergyDInternalEnergyRock(double phi,
                                                double phi0,
                                                double sl,
                                                double nl,
                                                double ul,
                                                double rho_r,
                                                double ur,
                                                double cv) const
{
  return cv * rho_r * (-phi0 + 1);
}

double
RichardsEnergyModel::DEnergyDCellVolume(double phi,
                                        double phi0,
                                        double sl,
                                        double nl,
                                        double ul,
                                        double rho_r,
                                        double ur,
                                        double cv) const
{
  return nl * phi * sl * ul + rho_r * ur * (-phi0 + 1);
}

} // namespace Relations
} // namespace Energy
} // namespace ATS_Physics
} // namespace Amanzi
