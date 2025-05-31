/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatsky (dasvyat@lanl.gov)
*/

//! Gravity driven preferential flow
/*!

Solves Preferential Flow equation:

.. math::
\frac{\partial \Theta}{\partial t} - \nabla \cdot ( k_m(s_m)^a \hat{z} + k \( \nabla p + \rho g \hat{z} ) ) = Q_w

.. _preferential-spec:
.. admonition:: preferential-spec

   * `"domain`" ``[string]`` **"domain"**  Defaults to the subsurface mesh.

   * `"primary variable key`" ``[string]`` The primary variable associated with
     this PK, typically `"DOMAIN-pressure`" Note there is no default -- this
     must be provided by the user.

   * `"boundary conditions`" ``[list]`` Defaults to Neuman, 0 normal
     flux.  See `Flow-specific Boundary Conditions`_

   * `"permeability type`" ``[string]`` **scalar** This controls the
     number of values needed to specify the absolute permeability.
     One of:

     - `"scalar`" Requires one scalar value.
     - `"horizontal and vertical`" Requires two values, horizontal
       then vertical.
     - `"diagonal tensor`" Requires dim values: {xx, yy} or {xx, yy,
       zz}
     - `"full tensor`". (Note symmetry is required.)  Either {xx, yy,
       xy} or {xx,yy,zz,xy,xz,yz}.

   * `"water retention evaluator`" evaluator-``[wrm-spec]`` The water
     retention curve.  This needs to go away, and should get moved to
     State.

   * `"water retention evaluator for gravity terms`" evaluator-``[wrm-spec]`` The water
     retention curve.  This needs to go away, and should get moved to
     State.


   IF
   * `"source term`" ``[bool]`` **false** Is there a source term?

   THEN
   * `"source key`" ``[string]`` **DOMAIN-water_source** Typically not
     set, as the default is good. ``[mol s^-1]``
   * `"source term is differentiable`" ``[bool]`` **true** Can the
     source term be differentiated with respect to the primary
     variable?
   * `"explicit source term`" ``[bool]`` **false** Apply the source
     term from the previous timestep.

   END

   Math and solver algorithm options:

   * `"diffusion`" ``[pde-diffusion-spec]`` The (forward) diffusion
     operator, see PDE_Diffusion_.

   * `"diffusion preconditioner`" ``[pde-diffusion-spec]``
     **optional** The inverse of the diffusion operator.  See
     PDE_Diffusion_.  Typically this is only needed to set Jacobian
     options, as all others probably should match those in
     `"diffusion`", and default to those values.

   * `"surface rel perm strategy`" ``[string]`` **none** Approach for
     specifying the relative permeabiilty on the surface face.
     `"clobber`" is frequently used for cases where a surface rel
     perm will be provided.  One of:

     - `"none`" : use the upwind direction to determine whether to
       use the boundary face or internal cell
     - `"clobber`" : always use the boundary face rel perm
     - `"max`" : use the max of the boundary face and internal cell
       values
     - `"unsaturated`" : Uses the boundary face when the internal
       cell is not saturated.

   * `"relative permeability method`" ``[string]`` **upwind with Darcy
     flux** Relative permeability is defined on cells, but must be
     calculated on faces to multiply a flux.  There are several
     methods commonly used.  Note these can significantly change
     answers -- you don't want to change these unless you know what
     they mean.  One of:

     - `"upwind with Darcy flux`" First-order upwind method that is
       most common
     - `"upwind with gravity`" Upwinds according to the gravitational
       flux direction
     - `"cell centered`" This corresponds to the harmonic mean, and is
        most accurate if the problem is always wet, but has issues
        when it is dry.
     - `"arithmetic mean`" Face value is the mean of the neighboring
        cells.  Not a good method.

   Globalization and other process-based hacks:

   * `"modify predictor with consistent faces`" ``[bool]`` **false** In a
     face+cell diffusion discretization, this modifies the predictor to make
     sure that faces, which are a DAE, are consistent with the predicted cells
     (i.e. face fluxes from each sides match).

   * `"modify predictor for flux BCs`" ``[bool]`` **false** Infiltration into
     dry ground can be hard on solvers -- this tries to do the local nonlinear
     problem to ensure that face pressures are consistent with the
     prescribed flux in a predictor.

   * `"modify predictor via water content`" ``[bool]`` **false** Modifies the
     predictor using the method of Krabbenhoft [??] paper.  Effectively does a
     change of variables, extrapolating not in pressure but in water content,
     then takes the smaller of the two extrapolants.

   * `"max valid change in saturation in a timestep [-]`" ``[double]`` **-1**
     Rejects timesteps whose max saturation change is greater than this value.
     This can be useful to ensure temporally resolved solutions.  Usually a
     good value is 0.1 or 0.2.

   * `"max valid change in ice saturation in a timestep [-]`" ``[double]``
     **-1** Rejects timesteps whose max ice saturation change is greater than
     this value.  This can be useful to ensure temporally resolved solutions.
     Usually a good value is 0.1 or 0.2.

   * `"limit correction to pressure change [Pa]`" ``[double]`` **-1** If > 0,
     this limits an iterate's max pressure change to this value.  Not usually
     helpful.

   * `"limit correction to pressure change when crossing atmospheric [Pa]`" ``[double]`` **-1**
     If > 0, this limits an iterate's max pressure change
     to this value when they cross atmospheric pressure.  Not usually helpful.

   Discretization / operators / solver controls:

   * `"accumulation preconditioner`" ``[pde-accumulation-spec]`` **optional**
     The inverse of the accumulation operator.  See PDE_Accumulation_.
     Typically not provided by users, as defaults are correct.

   * `"absolute error tolerance`" ``[double]`` **2750.0** ``[mol]``

   * `"compute boundary values`" ``[bool]`` **false** Used to include boundary
     face unknowns on discretizations that are cell-only (e.g. FV).  This can
     be useful for surface flow or other wierd boundary conditions.  Usually
     provided by MPCs that need them.

   Physics control:

   * `"permeability rescaling`" ``[double]`` **1e7** Typically 1e7 or order
     :math:`sqrt(K)` is about right.  This rescales things to stop from
     multiplying by small numbers (permeability) and then by large number
     (:math:`\rho / \mu`).

   IF
   * `"coupled to surface via flux`" ``[bool]`` **false** If true, apply
     surface boundary conditions from an exchange flux.  Note, if this is a
     coupled problem, it is probably set by the MPC.  No need for a user to
     set it.

   THEN
   * `"surface-subsurface flux key`" ``[string]`` **DOMAIN-surface_subsurface_flux**

   END

   * `"coupled to surface via head`" ``[bool]`` **false** If true, apply
     surface boundary conditions from the surface pressure (Dirichlet).

*/


/*
  Debugging including these..

   INCLUDES:
   - ``[pk-physical-bdf-default-spec]`` A `PK: Physical and BDF`_ spec.

   EVALUATORS:
   - `"conserved quantity`"
   - `"mass density`"
   - `"molar density`"
   - `"permeability`"
   - `"conductivity`"
   - `"saturation`"
   - `"primary variable`"


   Everything below this point is usually not provided by the user, but are
   documented here for completeness.

   KEYS:

   - `"conserved quantity`" **DOMAIN-water_content** Typically
      not set, as the default is good. ``[mol]``
   - `"mass density`" **DOMAIN-mass_density_liquid** liquid water
      density ``[kg m^-3]``
   - `"molar density`" **DOMAIN-molar_density_liquid** liquid
      water density ``[mol m^-3]``
   - `"permeability key`" **DOMAIN-permeability** permeability of the
      soil medium ``[m^2]``
   - `"conductivity key`" **DOMAIN-relative_permeability** scalar
      coefficient of the permeability ``[-]``
   - `"upwind conductivity key`" ``[string]``
      **DOMAIN-upwind_relative_permeability** upwinded (face-based) scalar
      coefficient of the permeability.  Note the units of this are strange, but
      this represents :math:`\frac{n_l k_r}{\mu}` ``[mol kg^-1 s^1 m^-2]``
   - `"darcy flux key`" **DOMAIN-water_flux** water flux across a face ``[mol s^-1]``
   - `"darcy flux direction key`" **DOMAIN-water_flux_direction**
      direction of the darcy flux (used in upwinding :math:`k_r`) ``[??]``
   - `"darcy velocity key`" **DOMAIN-darcy_velocity** darcy velocity
      vector, interpolated from faces to cells ``[m s^-1]``
   - `"saturation key`" **DOMAIN-saturation_liquid** volume
      fraction of the liquid phase ``[-]``
   - `"saturation gas key`" **DOMAIN-saturation_gas** volume
      fraction of the gas phase ``[-]``

*/


#ifndef PK_FLOW_PREFERENTIAL_HH_
#define PK_FLOW_PREFERENTIAL_HH_

#include "richards.hh"

#include "PDE_DiffusionFactory.hh"
#include "PDE_Accumulation.hh"
#include "PK_Factory.hh"
#include "pk_physical_bdf_default.hh"

namespace Amanzi {

// forward declarations
class MPCSubsurface;
class PredictorDelegateBCFlux;
class PrimaryVariableFieldEvaluator;
namespace WhetStone {
class Tensor;
}

namespace Flow {

class Preferential : public Richards {
 public:
  Preferential(Teuchos::ParameterList& pk_tree,
               const Teuchos::RCP<Teuchos::ParameterList>& plist,
               const Teuchos::RCP<State>& S,
               const Teuchos::RCP<TreeVector>& solution);

  // Virtual destructor
  virtual ~Preferential() {}

 protected:
  // Create of physical evaluators.
  virtual void SetupPhysicalEvaluators_() override;
  virtual void RequireNonlinearCoefficient_(const Key& key, const std::string& loc) override;

  // Builds tensor K, along with faced-based Krel if needed by the rel-perm method
  virtual bool UpdatePermeabilityData_(const Tag& tag) override;
  virtual bool UpdatePermeabilityDerivativeData_(const Tag& tag) override;


 protected:
  Key coef_grav_key_;
  Key dcoef_grav_key_;

 private:
  // factory registration
  static RegisteredPKFactory<Preferential> reg_;

  // Preferential has a friend in couplers...
  friend class Amanzi::MPCSubsurface;
};

} // namespace Flow
} // namespace Amanzi

#endif
