/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@ornl.gov)
*/

/*!

When including saltwater, the mass density of the liquid phase is a function of
the concentration (or mol ratio) of salt.  This is traditionally computed by a
known mass density at a reference value of the concentration of salt, e.g.

:math:`\rho_{ref} = \rho ( C_{ref} )`

In ATS, we compute transport in the primary variable of mol ratio [mol i / mol H2O],
so we convert units to write it as:

:math:`\rho( \xi^i ) = \rho_{f} + (\rho_{ref} - \rho_f) \frac{M^{salt} n}{C_{ref}) \xi^{salt}`

Note that the molar density provided [mols H2O / m^3] is a constant.

.. _eos-saltwater-spec
.. admonition:: eos-saltwater-spec

   * `"reference saltwater concentration [kg m^-3]`" ``[double]`` **35**
     Reference concentration of saltwater.

   * `"reference saltwater mass density [kg m^-3]`" ``[double]`` **1025** Mass
     density of saltwater at the reference concentration.

   * `"salt molar mass [g mol^-1]`" ``[double]`` **58.5** Mass of one mol of salt.

   * `"water molar mass [g mol^-1]`" ``[double]`` **18.0153** Mass of one mol of H2O.

   * `"fresh water mass density [kg m^-3]" ``[double]`` **1000** Mass density
     of fresh water.

*/

#ifndef AMANZI_RELATIONS_EOS_SW_HH_
#define AMANZI_RELATIONS_EOS_SW_HH_

#include "Teuchos_ParameterList.hpp"
#include "eos.hh"
#include "dbc.hh"

namespace Amanzi {
namespace Relations {

class EOS_SW : public EOS {
 public:
  explicit EOS_SW(Teuchos::ParameterList& eos_plist);
  virtual ~EOS_SW(){};

  // Virtual methods that form the EOS
  virtual double MassDensity(std::vector<double>& params) override;
  virtual double DMassDensityDMolarRatio(std::vector<double>& params) override;

  virtual double MolarDensity(std::vector<double>& params) override;
  virtual double DMolarDensityDMolarRatio(std::vector<double>& params) override;

  virtual bool IsTemperature() override { return false; }
  virtual bool IsPressure() override { return false; }
  virtual bool IsMolarRatio() override { return true; }

  // If molar mass is constant, we can take some shortcuts if we need both
  // molar and mass densities.  MolarMass() is undefined if
  // !IsConstantMolarMass()
  virtual bool IsConstantMolarMass() override { return false; }
  virtual double MolarMass() override
  {
    AMANZI_ASSERT(0);
    return 0.0;
  }

 protected:
  virtual void InitializeFromPlist_();

  Teuchos::ParameterList eos_plist_;
  double E_;
  double rho_f_;
  double n_l_;

 private:
  static Utils::RegisteredFactory<EOS, EOS_SW> factory_;
};

} // namespace Relations
} // namespace Amanzi

#endif
