/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Painter's permafrost model with freezing point depression, smoothed.
/*!

.. _wrm-fpd-smoothed-permafrost-spec:
.. admonition:: wrm-fpd-smoothed-permafrost-spec

    * `"reference temperature [K]`" ``[double]`` **273.15** The phase transition point
    * `"interfacial tension ice-water [mN m^-1]`" ``[double]`` **33.1**
    * `"interfacial tension air-water [mN m^-1]`" ``[double]`` **72.7**
    * `"smoothing width [K]`" ``[double]`` **1.** Smoothing out the freeze curve allows this to be slightly easier to solve.
    * `"latent heat [J kg^-1]`" ``[double]`` **3.34e5** Latent heat of fusion
    * `"water density [kg m^-3]`" ``[double]`` **998.87** Density of water.  Note this probably should use the calculated value.

 */

#ifndef AMANZI_FLOWRELATIONS_WRM_FPD_SMOOTHED_PERMAFROST_MODEL_
#define AMANZI_FLOWRELATIONS_WRM_FPD_SMOOTHED_PERMAFROST_MODEL_

#include "wrm_permafrost_model.hh"
#include "wrm_permafrost_factory.hh"

namespace Amanzi {
namespace Flow {

class WRM;

class WRMFPDSmoothedPermafrostModel : public WRMPermafrostModel {
 public:
  explicit WRMFPDSmoothedPermafrostModel(Teuchos::ParameterList& plist);

  // required methods from the base class
  // sats[0] = sg, sats[1] = sl, sats[2] = si
  virtual bool freezing(double T, double pc_liq, double pc_ice)
  {
    return pc_liq <= 0. ? T < 273.15 : pc_liq < pc_ice;
  }

  virtual void saturations(double pc_liq, double pc_ice, double (&sats)[3]);
  virtual void dsaturations_dpc_liq(double pc_liq, double pc_ice, double (&dsats)[3]);
  virtual void dsaturations_dpc_ice(double pc_liq, double pc_ice, double (&dsats)[3]);

 protected:
  double dp_;

 private:
  // factory registration
  static Utils::RegisteredFactory<WRMPermafrostModel, WRMFPDSmoothedPermafrostModel> factory_;
};


} // namespace Flow
} // namespace Amanzi

#endif
