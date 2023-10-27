/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Internal energy model for air and water vapor.
/*!

.. _iem-water-vapor-model-spec:
.. admonition:: iem-water-vapor-model-spec

    * `"latent heat [J mol^-1]`" ``[double]`` Latent heat of vaporization
    * `"heat capacity [J mol^-1 K^-1]`" ``[double]`` C_v

*/

#ifndef AMANZI_ENERGY_RELATIONS_IE_WATER_VAPOR_
#define AMANZI_ENERGY_RELATIONS_IE_WATER_VAPOR_

#include "Teuchos_ParameterList.hpp"

namespace Amanzi {
namespace Energy {

class IEMWaterVapor {
 public:
  IEMWaterVapor(Teuchos::ParameterList& plist);

  bool IsMolarBasis() { return true; }

  double InternalEnergy(double temp, double mol_frac_gas);
  double DInternalEnergyDT(double temp, double mol_frac_gas);
  double DInternalEnergyDomega(double temp, double mol_frac_gas);

 private:
  void InitializeFromPlist_();

  Teuchos::ParameterList plist_;

  double Cv_air_;            // units: MJ/(mol-K)
  double heat_vaporization_; // units: MJ/mol
};

} // namespace Energy
} // namespace Amanzi

#endif
