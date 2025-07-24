/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*

Interfrost model for saturated

 */

#ifndef AMANZI_FLOWRELATIONS_WRM_INTERFROST_PERMAFROST_MODEL_
#define AMANZI_FLOWRELATIONS_WRM_INTERFROST_PERMAFROST_MODEL_

#include "wrm_permafrost_model.hh"
#include "wrm_permafrost_factory.hh"

namespace Amanzi {
namespace Flow {

class WRM;

class WRMInterfrostPermafrostModel : public WRMPermafrostModel {
 public:
  explicit WRMInterfrostPermafrostModel(Teuchos::ParameterList& plist)
    : WRMPermafrostModel(plist)
  {
    W_ = plist.get<double>("W [K]");
    sr_ = plist.get<double>("residual saturation [-]");
    AMANZI_ASSERT(sr_ > 0.);
    AMANZI_ASSERT(sr_ < 1.);
  }

  // required methods from the base class
  // sats[0] = sg, sats[1] = sl, sats[2] = si
  virtual bool freezing(double T, double pc_liq, double pc_ice) { return T < 273.15; }

  virtual void saturations(double pc_liq, double pc_ice, double (&sats)[3]);
  virtual void dsaturations_dpc_liq(double pc_liq, double pc_ice, double (&dsats)[3]);
  virtual void dsaturations_dpc_ice(double pc_liq, double pc_ice, double (&dsats)[3]);

 protected:
  double W_;
  double sr_;

 private:
  // factory registration
  static Utils::RegisteredFactory<WRMPermafrostModel, WRMInterfrostPermafrostModel> factory_;
};


} // namespace Flow
} // namespace Amanzi

#endif
