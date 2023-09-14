/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Internal energy based on a quadratic fit to data.
/*!

Quadratic internal energy model -- function of Cv and temperature

.. math::

    u = L_f + C_v * (T - T_{ref}) + b(T - T_{ref})^2

`"IEM type`" = `"quadratic`"

.. _IEM-model-quadratic-spec:
.. admonition:: IEM-model-quadratic-spec

    * `"reference temperature [K]`" ``[double]`` **273.15** The phase transition point, T_ref above.

    ONE OF

    * `"latent heat [J kg^-1]`" ``[double]`` Latent heat of fusion, L_f above
    * `"heat capacity [J kg^-1 K^-1]`" ``[double]`` C_v above
    * `"quadratic b [J kg^-1 K^-2]`" ``[double]`` b above

    OR

    * `"latent heat [J mol^-1]`" ``[double]`` Latent heat of fusion, L_f above.
    * `"heat capacity [J mol^-1 K^-1]`" ``[double]`` C_v above
    * `"quadratic b [J mol^-1 K^-2]`" ``[double]`` b above

    END

*/

#ifndef AMANZI_ENERGYRELATIONS_IEM_QUADRATIC_
#define AMANZI_ENERGYRELATIONS_IEM_QUADRATIC_

#include "Teuchos_ParameterList.hpp"

#include "iem.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Energy {

class IEMQuadratic : public IEM {
 public:
  explicit IEMQuadratic(Teuchos::ParameterList& plist);

  bool IsMolarBasis() { return molar_basis_; }

  double InternalEnergy(double temp);
  double DInternalEnergyDT(double temp);

 private:
  virtual void InitializeFromPlist_();

  Teuchos::ParameterList plist_;

  double u0_;
  double ka_;
  double kb_;
  double T0_; // units: K
  bool molar_basis_;

 private:
  // iem factor registration
  static Utils::RegisteredFactory<IEM, IEMQuadratic> factory_;
};

} // namespace Energy
} // namespace Amanzi

#endif
