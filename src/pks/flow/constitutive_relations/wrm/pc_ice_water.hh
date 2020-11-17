/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! Capillary pressure of ice on water.

/*!
.. _pc-ice-water-spec
.. admonition:: pc-ice-water-spec

    * `"reference temperature [K]`" ``[double]`` **273.15** The phase transition point, T_ref above
    * `"interfacial tension ice-water [mN m^-1]`" ``[double]`` **33.1**
    * `"interfacial tension air-water [mN m^-1]`" ``[double]`` **72.7**
    * `"smoothing width [K]`" ``[double]`` **0.** Smoothing out the freeze curve allows this to be slightly easier to solve.

    ONE OF

    * `"latent heat [J kg^-1]`" ``[double]`` **3.34e5** Latent heat of fusion

    OR

    * `"latent heat [J mol^-1]`" ``[double]`` Latent heat of fusion

    END

*/

#ifndef AMANZI_FLOW_RELATIONS_PC_ICE_WATER_
#define AMANZI_FLOW_RELATIONS_PC_ICE_WATER_

#include "Teuchos_ParameterList.hpp"

namespace Amanzi {
namespace Flow {

class PCIceWater {

public:
  explicit PCIceWater(Teuchos::ParameterList& plist);

  // required methods from the base class
  bool IsMolarBasis() { return molar_basis_; }
  double CapillaryPressure(double T, double dens);
  double DCapillaryPressureDT(double T, double dens);
  double DCapillaryPressureDRho(double T, double dens);

private:
  void InitializeFromPlist_();

  Teuchos::ParameterList pc_plist_;

  bool molar_basis_;
  double sigma_ice_liq_;
  double sigma_gas_liq_;
  double heat_fusion_;
  double gamma_;
  double T0_;
  double halfwidth_;
};

} //namespace
} //namespace

#endif
