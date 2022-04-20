/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! An advection-diffusion equation for surface energy in two phases.

/*!

This is simply a surface energy equation that places a few more requirements
on the base class.  It could probably go away if we refactor to remove
hard-coded evaluators.

.. _energy-surface-ice-pk-spec:
.. admonition:: energy-surface-ice-pk-spec

    These are typically not set by the user:

    * `"coupled to subsurface via temperature`" ``[bool]`` **false** A coupling
      scheme, provided by MPC.

    * `"coupled to subsurface via flux`" ``[bool]`` **false** A coupling
      scheme, provided by MPC.

    * `"subsurface domain name`" ``[string]`` **optional** If one of the above
      coupling schemes is turned on, we need to know the subsurface mesh.
      Provided by MPC.

    INCLUDES:

    - ``[energy-pk-spec]``  See `Energy Base PK`_

*/


#ifndef PKS_ENERGY_SURFACE_ICE_HH_
#define PKS_ENERGY_SURFACE_ICE_HH_

#include "PK_Factory.hh"
#include "energy_base.hh"

namespace Amanzi {

// forward declarations
class Function;

namespace Energy {

class EnergySurfaceIce : public EnergyBase {
 public:
  EnergySurfaceIce(Teuchos::ParameterList& FElist,
                   const Teuchos::RCP<Teuchos::ParameterList>& plist,
                   const Teuchos::RCP<State>& S,
                   const Teuchos::RCP<TreeVector>& solution);

  // -- Initialize owned (dependent) variables.
  virtual void Initialize(const Teuchos::Ptr<State>& S);

 protected:
  // -- setup the evaluators
  virtual void SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S);

  // -- get enthalpy as a function of Dirichlet boundary data.  Note that this
  //    will get replaced by a better system when we get maps on the boundary
  //    faces.
  //  virtual void ApplyDirichletBCsToEnthalpy_(const Teuchos::Ptr<State>& S);

  virtual void AddSources_(const Teuchos::Ptr<State>& S,
                           const Teuchos::Ptr<CompositeVector>& f);
  virtual void AddSourcesToPrecon_(const Teuchos::Ptr<State>& S, double h);

 protected:
  // simple heat condution term, q = K_s2a * (Tair - Tsurf)
  // air temperature function of time (not space)
  double K_surface_to_air_;

  // flags
  bool standalone_mode_;
  bool is_energy_source_term_;
  bool is_water_source_term_;
  bool is_air_conductivity_;
  bool water_source_only_if_unfrozen_;

  // keys
  Key domain_ss_;

 private:
  // factory registration
  static RegisteredPKFactory<EnergySurfaceIce> reg_;

};

} // namespace Energy
} // namespace Amanzi

#endif
