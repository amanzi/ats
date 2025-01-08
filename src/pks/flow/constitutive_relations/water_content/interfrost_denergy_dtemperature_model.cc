/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  The interfrost denergy_dtemperature model is an algebraic model with dependencies.

  Generated via evaluator_generator with:
Interfrost water content portion sl.

*/

#include "Teuchos_ParameterList.hpp"
#include "dbc.hh"
#include "interfrost_denergy_dtemperature_model.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

// Constructor from ParameterList
InterfrostDenergyDtemperatureModel::InterfrostDenergyDtemperatureModel(
  Teuchos::ParameterList& plist)
{
  InitializeFromPlist_(plist);
}


// Initialize parameters
void
InterfrostDenergyDtemperatureModel::InitializeFromPlist_(Teuchos::ParameterList& plist)
{
  W_ = plist.get<double>("W [K]");
}


// main method
double
InterfrostDenergyDtemperatureModel::DEnergyDTCoef(double phi,
                                                  double sl,
                                                  double nl,
                                                  double si,
                                                  double ni,
                                                  double rhos,
                                                  double T) const
{
  return 0.0060171102 * ni * phi *
           ((T >= 273.15) ?
              (0.0) :
              (-0.95 * (2 * T - 546.3) * exp(-pow(T - 273.15, 2) / pow(W_, 2)) / pow(W_, 2))) +
         1.0e-6 * phi * (37.111518 * ni * si + 75.3399846 * nl * sl) + 0.000835 * rhos * (-phi + 1);
}

double
InterfrostDenergyDtemperatureModel::DDEnergyDTCoefDPorosity(double phi,
                                                            double sl,
                                                            double nl,
                                                            double si,
                                                            double ni,
                                                            double rhos,
                                                            double T) const
{
  return 3.7111518e-5 * ni * si +
         0.0060171102 * ni *
           ((T >= 273.15) ?
              (0.0) :
              (-0.95 * (2 * T - 546.3) * exp(-pow(T - 273.15, 2) / pow(W_, 2)) / pow(W_, 2))) +
         7.53399846e-5 * nl * sl - 0.000835 * rhos;
}

double
InterfrostDenergyDtemperatureModel::DDEnergyDTCoefDSaturationLiquid(double phi,
                                                                    double sl,
                                                                    double nl,
                                                                    double si,
                                                                    double ni,
                                                                    double rhos,
                                                                    double T) const
{
  return 7.53399846e-5 * nl * phi;
}

double
InterfrostDenergyDtemperatureModel::DDEnergyDTCoefDMolarDensityLiquid(double phi,
                                                                      double sl,
                                                                      double nl,
                                                                      double si,
                                                                      double ni,
                                                                      double rhos,
                                                                      double T) const
{
  return 7.53399846e-5 * phi * sl;
}

double
InterfrostDenergyDtemperatureModel::DDEnergyDTCoefDSaturationIce(double phi,
                                                                 double sl,
                                                                 double nl,
                                                                 double si,
                                                                 double ni,
                                                                 double rhos,
                                                                 double T) const
{
  return 3.7111518e-5 * ni * phi;
}

double
InterfrostDenergyDtemperatureModel::DDEnergyDTCoefDMolarDensityIce(double phi,
                                                                   double sl,
                                                                   double nl,
                                                                   double si,
                                                                   double ni,
                                                                   double rhos,
                                                                   double T) const
{
  return 3.7111518e-5 * phi * si +
         0.0060171102 * phi *
           ((T >= 273.15) ?
              (0.0) :
              (-0.95 * (2 * T - 546.3) * exp(-pow(T - 273.15, 2) / pow(W_, 2)) / pow(W_, 2)));
}

double
InterfrostDenergyDtemperatureModel::DDEnergyDTCoefDDensityRock(double phi,
                                                               double sl,
                                                               double nl,
                                                               double si,
                                                               double ni,
                                                               double rhos,
                                                               double T) const
{
  return -0.000835 * phi + 0.000835;
}

double
InterfrostDenergyDtemperatureModel::DDEnergyDTCoefDTemperature(double phi,
                                                               double sl,
                                                               double nl,
                                                               double si,
                                                               double ni,
                                                               double rhos,
                                                               double T) const
{
  return 0.0060171102 * ni * phi *
         ((T >= 273.15) ?
            (0) :
            (-1.9 * exp(-pow(T - 273.15, 2) / pow(W_, 2)) / pow(W_, 2) +
             0.95 * pow(2 * T - 546.3, 2) * exp(-pow(T - 273.15, 2) / pow(W_, 2)) / pow(W_, 4)));
}

} // namespace Relations
} // namespace Flow
} // namespace Amanzi
