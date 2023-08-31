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

#ifndef AMANZI_ENERGY_LIQUID_ICE_ENERGY_MODEL_HH_
#define AMANZI_ENERGY_LIQUID_ICE_ENERGY_MODEL_HH_

namespace Amanzi {
namespace Energy {
namespace Relations {

class LiquidIceEnergyModel {
 public:
  explicit LiquidIceEnergyModel(Teuchos::ParameterList& plist);

  double Energy(double phi,
                double phi0,
                double sl,
                double nl,
                double ul,
                double si,
                double ni,
                double ui,
                double rho_r,
                double ur,
                double cv) const;

  double DEnergyDPorosity(double phi,
                          double phi0,
                          double sl,
                          double nl,
                          double ul,
                          double si,
                          double ni,
                          double ui,
                          double rho_r,
                          double ur,
                          double cv) const;
  double DEnergyDBasePorosity(double phi,
                              double phi0,
                              double sl,
                              double nl,
                              double ul,
                              double si,
                              double ni,
                              double ui,
                              double rho_r,
                              double ur,
                              double cv) const;
  double DEnergyDSaturationLiquid(double phi,
                                  double phi0,
                                  double sl,
                                  double nl,
                                  double ul,
                                  double si,
                                  double ni,
                                  double ui,
                                  double rho_r,
                                  double ur,
                                  double cv) const;
  double DEnergyDMolarDensityLiquid(double phi,
                                    double phi0,
                                    double sl,
                                    double nl,
                                    double ul,
                                    double si,
                                    double ni,
                                    double ui,
                                    double rho_r,
                                    double ur,
                                    double cv) const;
  double DEnergyDInternalEnergyLiquid(double phi,
                                      double phi0,
                                      double sl,
                                      double nl,
                                      double ul,
                                      double si,
                                      double ni,
                                      double ui,
                                      double rho_r,
                                      double ur,
                                      double cv) const;
  double DEnergyDSaturationIce(double phi,
                               double phi0,
                               double sl,
                               double nl,
                               double ul,
                               double si,
                               double ni,
                               double ui,
                               double rho_r,
                               double ur,
                               double cv) const;
  double DEnergyDMolarDensityIce(double phi,
                                 double phi0,
                                 double sl,
                                 double nl,
                                 double ul,
                                 double si,
                                 double ni,
                                 double ui,
                                 double rho_r,
                                 double ur,
                                 double cv) const;
  double DEnergyDInternalEnergyIce(double phi,
                                   double phi0,
                                   double sl,
                                   double nl,
                                   double ul,
                                   double si,
                                   double ni,
                                   double ui,
                                   double rho_r,
                                   double ur,
                                   double cv) const;
  double DEnergyDDensityRock(double phi,
                             double phi0,
                             double sl,
                             double nl,
                             double ul,
                             double si,
                             double ni,
                             double ui,
                             double rho_r,
                             double ur,
                             double cv) const;
  double DEnergyDInternalEnergyRock(double phi,
                                    double phi0,
                                    double sl,
                                    double nl,
                                    double ul,
                                    double si,
                                    double ni,
                                    double ui,
                                    double rho_r,
                                    double ur,
                                    double cv) const;
  double DEnergyDCellVolume(double phi,
                            double phi0,
                            double sl,
                            double nl,
                            double ul,
                            double si,
                            double ni,
                            double ui,
                            double rho_r,
                            double ur,
                            double cv) const;

 protected:
  void InitializeFromPlist_(Teuchos::ParameterList& plist);

 protected:
};

} // namespace Relations
} // namespace Energy
} // namespace Amanzi

#endif
