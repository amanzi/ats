/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/
/*!

Painter's permafrost model with freezing point depression.

`"permafrost wrm type`" = `"fpd permafrost model`"

.. _permafrost-wrm-fpd-permafrost-model-spec:
.. admonition:: permafrost-wrm-fpd-permafrost-model-spec

    * `"empty`" ``[double]`` **-1** No parameters are required for this model

*/

#ifndef AMANZI_FLOWRELATIONS_WRM_FPD_PERMAFROST_MODEL_
#define AMANZI_FLOWRELATIONS_WRM_FPD_PERMAFROST_MODEL_

#include "wrm_permafrost_model.hh"
#include "wrm_permafrost_factory.hh"

namespace Amanzi {
namespace Flow {

class WRM;

class WRMFPDPermafrostModel : public WRMPermafrostModel {
 public:
  explicit WRMFPDPermafrostModel(Teuchos::ParameterList& plist) : WRMPermafrostModel(plist) {}

  // required methods from the base class
  // sats[0] = sg, sats[1] = sl, sats[2] = si
  virtual bool freezing(double T, double pc_liq, double pc_ice)
  {
    return pc_liq <= 0. ? pc_ice > 0. : pc_liq < pc_ice;
  }

  virtual void saturations(double pc_liq, double pc_ice, double (&sats)[3]);
  virtual void dsaturations_dpc_liq(double pc_liq, double pc_ice, double (&dsats)[3]);
  virtual void dsaturations_dpc_ice(double pc_liq, double pc_ice, double (&dsats)[3]);

 protected:
  double deriv_regularization_;

 private:
  // factory registration
  static Utils::RegisteredFactory<WRMPermafrostModel, WRMFPDPermafrostModel> factory_;
};


} // namespace Flow
} // namespace Amanzi

#endif
