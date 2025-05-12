/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*!

Internal energy model for air and water vapor.

.. math::

   u = (1 + 0.622 \omega) C_v^{air} (T - 273.15) + \omega L_v

.. _iem-water-vapor-spec:
.. admonition:: iem-water-vapor-spec

   * `"latent heat [J mol^-1]`" ``[double]`` Latent heat of vaporization,
     :math:`L_v`
   * `"heat capacity [J mol^-1 K^-1]`" ``[double]`` Heat capacity of air,
     :math:`C_v^{air}`

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
