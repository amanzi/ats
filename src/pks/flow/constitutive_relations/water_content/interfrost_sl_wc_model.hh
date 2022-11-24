/*
  The interfrost sl water content model is an algebraic model with dependencies.

  Generated via evaluator_generator with:
Interfrost water content portion sl.
    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOW_INTERFROST_SL_WC_MODEL_HH_
#define AMANZI_FLOW_INTERFROST_SL_WC_MODEL_HH_

namespace Amanzi {
namespace Flow {
namespace Relations {

class InterfrostSlWcModel {
 public:
  explicit InterfrostSlWcModel(Teuchos::ParameterList& plist);

  double WaterContent(double phi, double sl, double nl, double ni, double cv) const;

  double DWaterContentDPorosity(double phi, double sl, double nl, double ni, double cv) const;
  double
  DWaterContentDSaturationLiquid(double phi, double sl, double nl, double ni, double cv) const;
  double
  DWaterContentDMolarDensityLiquid(double phi, double sl, double nl, double ni, double cv) const;
  double
  DWaterContentDMolarDensityIce(double phi, double sl, double nl, double ni, double cv) const;
  double DWaterContentDCellVolume(double phi, double sl, double nl, double ni, double cv) const;

 protected:
  void InitializeFromPlist_(Teuchos::ParameterList& plist);

 protected:
};

} // namespace Relations
} // namespace Flow
} // namespace Amanzi

#endif
