/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*!

This model is based on the emperical freezing curve used by Sutra-Ice,
documented in papers by Voss & Walvoord.

.. _permafrost-wrm-sutra-permafrost-model-spec:
.. admonition:: permafrost-wrm-sutra-permafrost-model-spec

   * `"temperature transition [K]`" ``[double]`` thickness of the transition from frozen to thawed
   * `"residual saturation [-]`" ``[double]`` Standard residual saturation
   * `"freezing point [K]`" ``[double]`` **273.15**


 */

#ifndef AMANZI_FLOWRELATIONS_WRM_SUTRA_PERMAFROST_MODEL_
#define AMANZI_FLOWRELATIONS_WRM_SUTRA_PERMAFROST_MODEL_

#include "wrm_permafrost_model.hh"
#include "wrm_permafrost_factory.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {

class WRM;

class WRMSutraPermafrostModel : public WRMPermafrostModel {
 public:
  explicit WRMSutraPermafrostModel(Teuchos::ParameterList& plist)
    : WRMPermafrostModel(plist)
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
} // namespace ATS_Physics
} // namespace Amanzi

#endif
