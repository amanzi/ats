/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*

Painter's permafrost model.

 */

#ifndef AMANZI_FLOWRELATIONS_WRM_OLD_PERMAFROST_MODEL_
#define AMANZI_FLOWRELATIONS_WRM_OLD_PERMAFROST_MODEL_

#include "wrm_permafrost_model.hh"
#include "wrm_permafrost_factory.hh"

namespace Amanzi {
namespace Flow {

class WRM;

class WRMOldPermafrostModel : public WRMPermafrostModel {
 public:
  WRMOldPermafrostModel(Teuchos::ParameterList& plist)
    : WRMPermafrostModel(plist)
  {}

  // required methods from the base class
  // sats[0] = s_g, sats[1] = s_l, sats[2] = s_i
  virtual bool freezing(double T, double pc_liq, double pc_ice) { return T < 273.15; }
  virtual void saturations(double pc_liq, double pc_ice, double (&sats)[3]);
  virtual void dsaturations_dpc_liq(double pc_liq, double pc_ice, double (&dsats)[3]);
  virtual void dsaturations_dpc_ice(double pc_liq, double pc_ice, double (&dsats)[3]);

 private:
  static Utils::RegisteredFactory<WRMPermafrostModel, WRMOldPermafrostModel> factory_;
};


} // namespace Flow
} // namespace Amanzi

#endif
