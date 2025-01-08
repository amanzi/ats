/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  The surface ice energy model is an algebraic model with dependencies.

  Generated via evaluator_generator with:
Energy evaulator for ice+liquid surface water.

*/

#ifndef AMANZI_ENERGY_SURFACE_ICE_ENERGY_MODEL_HH_
#define AMANZI_ENERGY_SURFACE_ICE_ENERGY_MODEL_HH_

namespace Amanzi {
namespace Energy {
namespace Relations {

class SurfaceIceEnergyModel {
 public:
  explicit SurfaceIceEnergyModel(Teuchos::ParameterList& plist);

  double Energy(double h, double eta, double nl, double ul, double ni, double ui, double cv) const;

  double
  DEnergyDPondedDepth(double h, double eta, double nl, double ul, double ni, double ui, double cv)
    const;
  double DEnergyDUnfrozenFraction(double h,
                                  double eta,
                                  double nl,
                                  double ul,
                                  double ni,
                                  double ui,
                                  double cv) const;
  double DEnergyDMolarDensityLiquid(double h,
                                    double eta,
                                    double nl,
                                    double ul,
                                    double ni,
                                    double ui,
                                    double cv) const;
  double DEnergyDInternalEnergyLiquid(double h,
                                      double eta,
                                      double nl,
                                      double ul,
                                      double ni,
                                      double ui,
                                      double cv) const;
  double DEnergyDMolarDensityIce(double h,
                                 double eta,
                                 double nl,
                                 double ul,
                                 double ni,
                                 double ui,
                                 double cv) const;
  double DEnergyDInternalEnergyIce(double h,
                                   double eta,
                                   double nl,
                                   double ul,
                                   double ni,
                                   double ui,
                                   double cv) const;
  double
  DEnergyDCellVolume(double h, double eta, double nl, double ul, double ni, double ui, double cv)
    const;

 protected:
  void InitializeFromPlist_(Teuchos::ParameterList& plist);

 protected:
};

} // namespace Relations
} // namespace Energy
} // namespace Amanzi

#endif
