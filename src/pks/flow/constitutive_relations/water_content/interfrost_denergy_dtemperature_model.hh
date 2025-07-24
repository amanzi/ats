/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  The interfrost denergy_dtemperature model is an algebraic model with dependencies.

  Generated via evaluator_generator with:
Interfrost water content portion sl.

*/

#ifndef AMANZI_FLOW_INTERFROST_DENERGY_DTEMPERATURE_MODEL_HH_
#define AMANZI_FLOW_INTERFROST_DENERGY_DTEMPERATURE_MODEL_HH_

namespace Amanzi {
namespace Flow {
namespace Relations {

class InterfrostDenergyDtemperatureModel {
 public:
  explicit InterfrostDenergyDtemperatureModel(Teuchos::ParameterList& plist);

  double DEnergyDTCoef(double phi,
                       double sl,
                       double nl,
                       double si,
                       double ni,
                       double rhos,
                       double T) const;

  double DDEnergyDTCoefDPorosity(double phi,
                                 double sl,
                                 double nl,
                                 double si,
                                 double ni,
                                 double rhos,
                                 double T) const;
  double DDEnergyDTCoefDSaturationLiquid(double phi,
                                         double sl,
                                         double nl,
                                         double si,
                                         double ni,
                                         double rhos,
                                         double T) const;
  double DDEnergyDTCoefDMolarDensityLiquid(double phi,
                                           double sl,
                                           double nl,
                                           double si,
                                           double ni,
                                           double rhos,
                                           double T) const;
  double DDEnergyDTCoefDSaturationIce(double phi,
                                      double sl,
                                      double nl,
                                      double si,
                                      double ni,
                                      double rhos,
                                      double T) const;
  double DDEnergyDTCoefDMolarDensityIce(double phi,
                                        double sl,
                                        double nl,
                                        double si,
                                        double ni,
                                        double rhos,
                                        double T) const;
  double DDEnergyDTCoefDDensityRock(double phi,
                                    double sl,
                                    double nl,
                                    double si,
                                    double ni,
                                    double rhos,
                                    double T) const;
  double DDEnergyDTCoefDTemperature(double phi,
                                    double sl,
                                    double nl,
                                    double si,
                                    double ni,
                                    double rhos,
                                    double T) const;

 protected:
  void InitializeFromPlist_(Teuchos::ParameterList& plist);

 protected:
  double W_;
};

} // namespace Relations
} // namespace Flow
} // namespace Amanzi

#endif
