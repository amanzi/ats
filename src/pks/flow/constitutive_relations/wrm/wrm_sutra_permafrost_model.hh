/*
Author: Ethan Coon

Sutra model for saturated

 */

#ifndef AMANZI_FLOWRELATIONS_WRM_SUTRA_PERMAFROST_MODEL_
#define AMANZI_FLOWRELATIONS_WRM_SUTRA_PERMAFROST_MODEL_

#include "wrm_permafrost_model.hh"
#include "wrm_permafrost_factory.hh"

namespace Amanzi {
namespace Flow {

class WRM;

class WRMSutraPermafrostModel : public WRMPermafrostModel {
 public:
  explicit WRMSutraPermafrostModel(Teuchos::ParameterList& plist) : WRMPermafrostModel(plist)
  {
    T0_ = plist.get<double>("freezing point [K]", 273.15);
    dT_ = plist.get<double>("temperature transition [K]");
    AMANZI_ASSERT(dT_ >= 0.);
    sr_ = plist.get<double>("residual saturation [-]");
    AMANZI_ASSERT(sr_ > 0.);
    AMANZI_ASSERT(sr_ < 1.);
  }

  // required methods from the base class
  // sats[0] = sg, sats[1] = sl, sats[2] = si
  virtual bool freezing(double T, double pc_liq, double pc_ice) { return T < T0_; }

  virtual void saturations(double pc_liq, double pc_ice, double (&sats)[3]);
  virtual void dsaturations_dpc_liq(double pc_liq, double pc_ice, double (&dsats)[3]);
  virtual void dsaturations_dpc_ice(double pc_liq, double pc_ice, double (&dsats)[3]);

 protected:
  double dT_;
  double T0_;
  double sr_;

 private:
  // factory registration
  static Utils::RegisteredFactory<WRMPermafrostModel, WRMSutraPermafrostModel> factory_;
};


} // namespace Flow
} // namespace Amanzi

#endif
