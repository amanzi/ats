/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatsky (dasvyat@lanl.gov)
*/

/*!

This PK solves sediment transport in surface flows. The assumpotion is that there is only one sediment component which
is called 'sediment'. The advection-diffusion equation for this component

.. math::
  \frac{\partial (\Theta \chi_i)}{\partial t} = - \boldsymbol{\nabla} \cdot (\boldsymbol{q} \chi_i)
  + \boldsymbol{\nabla} \cdot ( \tau \, \boldsymbol{\nabla} \chi_i) + Q_e - Q_t - Q_s,


As well as in transport ATS PK the primary variable is the **mole fraction** of a sediment with units of
[mol i mol liquid^-1]

   * `"sediment diffusion coefficient [m^2 s^-1]`" ``[molecular-diffusivity-spec]`` See
     below.

   Math and solver algorithm options:

   * `"diffusion`" ``[pde-diffusion-typedinline-spec]`` Diffusion drives the
     distribution.  Typically we use finite volume here, but mimetic schemes
     may be used.  See :ref:`Diffusion`


   * `"diffusion preconditioner`" ``[pde-diffusion-typedinline-spec]`` Inverse
     of the above.  Likely only Jacobian term options are needed here, as the
     others default to the same as the `"diffusion`" list.  See
     :ref:`Diffusion`.

   * `"inverse`" ``[inverse-typed-spec]`` :ref:`Inverse` method for the
     diffusion-dispersion solve.  See :ref:`Inverse`.

   * `"cfl`" ``[double]`` **1.0** Courant-Friedrichs-Lewy condition, a limiter
     on the timestep size relative to spatial size.  Must be <= 1.

   * `"advection spatial discretization order`" ``[int]`` **1** Defines
     accuracy of the spatial discretization in the advection term.  It permits
     values 1 or 2. Default value is 1 (donor upwind) but 2 (a limiter scheme)
     is much less numerically diffusive, and recommended for most cases.

   * `"temporal discretization order`" ``[int]`` **1** Defines accuracy of
     temporal discretization.  It permits values 1 (forward Euler) and 2 (a
     Runga-Kutta scheme).

   * `"reconstruction`" ``[reconstruction-spec]`` collects reconstruction
     parameters for use in reconstructing the velocity field for 2nd order
     advection schemes.  See :ref:`Reconstructions`.

   KEYS
   - `"primary variable`" **"mole_fraction"** [mol C mol H2O^-1]
   - `"liquid water content`" **"water_content"** This variable is a multiplier
     in in the accumulation term. This is often just `"water_content`", but
     not in frozen cases, where it must be defined by the user (typically as
     :math:`\phi s n_l |V|` in the subsurface, or :math:`h n_l |V|` on the
     surface.
   - `"molar density liquid`" [mol H2O m^-3]
   - `"water flux`" The face-based water flux in [mol H2O s^-1].
   - `"water source`" Defines the water injection rate [mol H2O m^-2 s^-1] in
     surface and [mol H2O m^-3 s^-1] in subsurface) which applies to
     concentrations specified by the `"geochemical conditions`".  Note that if
     this PK is coupled to a surface flow PK, the unit of the water source
     there *must* be in [mol m^-2 s^-1], *not* in [m s^-1] as is an option for
     that PK (e.g. `"water source in meters`" must be set to `"false`" in the
     overland flow PK).

     The sources for sediments transport are defined via "erosion rate", :math:`Q_e`, "traping rate", :math:`Q_t`, "settling rate", :math:`Q_s`
     evaluators. These evaluators are defined in unit of [m/s]. Moreover sediment density has to be defined as a scalar value.

*/

#ifndef AMANZI_ATS_SEDIMENTTRANSPORT_PK_HH_
#define AMANZI_ATS_SEDIMENTTRANSPORT_PK_HH_

#include <string>

// Amanzi
#include "PK_Factory.hh"
#include "State.hh"

// Transport
#include "transport_ats.hh"

/* ******************************************************************
The transport PK receives a reduced (optional) copy of a physical
state at time n and returns a different state at time n+1.

Unmodified physical quantaties in the returned state are the smart
pointers to the original variables.
****************************************************************** */

namespace Amanzi {
namespace Transport {


class SedimentTransport_PK : public Transport_ATS {
 public:
  SedimentTransport_PK(Teuchos::ParameterList& pk_tree,
                       const Teuchos::RCP<Teuchos::ParameterList>& glist,
                       const Teuchos::RCP<State>& S,
                       const Teuchos::RCP<TreeVector>& soln);

  void parseParameterList() override;
  void Initialize() override;

 protected:
  void SetupPhysicalEvaluators_() override;
  void AddSourceTerms_(double t0, double t1, Epetra_MultiVector& conserve_qty) override;
  ;


  Key sd_trapping_key_, sd_settling_key_, sd_erosion_key_, horiz_mixing_key_, sd_organic_key_;
  Key elevation_increase_key_;
  Key porosity_key_;
  Key plant_area_key_, stem_diameter_key_, stem_height_key_, stem_density_key_;

 private:
  double sediment_density_;


 private:
  // factory registration
  static RegisteredPKFactory<SedimentTransport_PK> reg_;
};

} // namespace Transport
} // namespace Amanzi

#endif
