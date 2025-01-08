/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  The richards water content model is an algebraic model with dependencies.

  Generated via evaluator_generator with:
Richards water content evaluator: the standard form as a function of liquid saturation.

*/

#ifndef AMANZI_FLOW_RICHARDS_WATER_CONTENT_MODEL_HH_
#define AMANZI_FLOW_RICHARDS_WATER_CONTENT_MODEL_HH_

namespace Amanzi {
namespace Flow {
namespace Relations {

class RichardsWaterContentModel {
 public:
  explicit RichardsWaterContentModel(Teuchos::ParameterList& plist);

  double WaterContent(double phi, double sl, double nl, double cv) const;

  double DWaterContentDPorosity(double phi, double sl, double nl, double cv) const;
  double DWaterContentDSaturationLiquid(double phi, double sl, double nl, double cv) const;
  double DWaterContentDMolarDensityLiquid(double phi, double sl, double nl, double cv) const;
  double DWaterContentDCellVolume(double phi, double sl, double nl, double cv) const;

 protected:
  void InitializeFromPlist_(Teuchos::ParameterList& plist);

 protected:
};

} // namespace Relations
} // namespace Flow
} // namespace Amanzi

#endif
