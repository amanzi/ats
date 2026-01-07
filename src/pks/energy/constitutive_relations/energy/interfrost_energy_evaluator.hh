/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Energy content evaluator for a liquid, ice system including the surrounding soil considering liquid compressibility.
/*!

Calculates energy, in [KJ], via the equation:

.. math::
  E = V * ( \phi (u_l s_l n_l (1 + \beta (pressure - 101325) ) + u_i s_i n_i)  + (1-\phi_0) u_r \rho_r )

`"evaluator type`" = `"interfrost energy`"

Note this equation assumes that porosity is compressible, but is based on the
uncompressed rock grain density (not bulk density).  This means that porosity
is the base, uncompressible value when used with the energy in the grain, but
the larger, compressible value when used with the energy in the water.

.. _evaluator-interfrost-energy-spec:
.. admonition:: evaluator-interfrost-energy-spec

   DEPENDENCIES:

   - `"porosity`" The porosity, including any compressibility. [-]
   - `"base porosity`" The uncompressed porosity (note this may be the same as
     porosity for incompressible cases) [-]
   - `"molar density liquid`" [mol m^-3]
   - `"saturation liquid`" [-]
   - `"internal energy liquid`" [KJ mol^-1]
   - `"molar density ice`" [mol m^-3]
   - `"saturation ice`" [-]
   - `"internal energy ice`" [KJ mol^-1]
   - `"density rock`" Units may be either [kg m^-3] or [mol m^-3]
   - `"internal energy rock`" Units may be either [KJ kg^-1] or [KJ mol^-1],
     but must be consistent with the above density.
   - `"cell volume`" [m^3]
   - `"pressure`" [Pa]

*/


#ifndef AMANZI_INTERFROST_ENERGY_EVALUATOR_HH_
#define AMANZI_INTERFROST_ENERGY_EVALUATOR_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Energy {
namespace Relations {

class InterfrostEnergyEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit InterfrostEnergyEvaluator(Teuchos::ParameterList& energy_plist);

  virtual Teuchos::RCP<Evaluator> Clone() const override;

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

 protected:
  Key phi_key_;
  Key phi0_key_;
  Key sl_key_;
  Key nl_key_;
  Key ul_key_;
  Key si_key_;
  Key ni_key_;
  Key ui_key_;
  Key rho_r_key_;
  Key ur_key_;
  Key cv_key_;
  Key pres_key_;

  double beta_;

 private:
  static Utils::RegisteredFactory<Evaluator, InterfrostEnergyEvaluator> reg_;
};

} // namespace Relations
} // namespace Energy
} // namespace Amanzi

#endif
