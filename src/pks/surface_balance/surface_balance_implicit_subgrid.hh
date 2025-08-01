/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
/*!

This is a balance PK whose conserved quantity is snow water equivalent.  The
energy balance comes in as it provides the energy needed to melt snow.  So
source terms include snow precipitation and snowmelt.  It also manages snow
density, which should get rethought a bit.

There is also some wierd hackiness here about area fractions -- see ATS Issue
#8.

`"PK type`" = `"surface balance implicit subgrid`"

.. _pk-surface-balance-implicit-subgrid-spec:
.. admonition:: pk-surface-balance-implicit-subgrid-spec

   * `"absolute error tolerance`" ``[double]`` **0.01** ``[m]`` A small amount
     of snow.

   INCLUDES:

   - ``[balance-pk-spec]`` This *is a* `Balance Equation`_

   KEYS:

   - `"conserved quantity key`" **DOMAIN-snow_water_equivalent** Sets the
     default conserved quantity key, so this is likely not supplied by the
     user. `[m]`
   - `"snow density key`" **DOMAIN-density** Default snow density key. `[kg
     m^-3]`
   - `"snow age key`" **DOMAIN-age** Default snow age key. `[d]`
   - `"new snow key`" **DOMAIN-source** Default new snow key. `[m SWE s^-1]`
   - `"area fractions key`" **DOMAIN-fractional_areas** Subgrid model
     fractional areas, see note above. `[-]`
   - `"snow death rate key`" **DOMAIN-death_rate** Deals with last tiny bit of
     snowmelt.

*/


#ifndef PK_SURFACE_BALANCE_IMPLICIT_SUBGRID_HH_
#define PK_SURFACE_BALANCE_IMPLICIT_SUBGRID_HH_

#include "PK_Factory.hh"
#include "surface_balance_base.hh"


namespace Amanzi {
namespace SurfaceBalance {

class ImplicitSubgrid : public SurfaceBalanceBase {
 public:
  ImplicitSubgrid(Teuchos::ParameterList& pk_tree,
                  const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                  const Teuchos::RCP<State>& S,
                  const Teuchos::RCP<TreeVector>& solution);

  // call to allow a PK to modify its own list or lists of its children.
  virtual void parseParameterList() override;

  // main methods
  // -- Setup data.
  virtual void Setup() override;

  // -- Initialize owned (dependent) variables.
  virtual void Initialize() override;

  virtual bool ModifyPredictor(double h,
                               Teuchos::RCP<const TreeVector> u0,
                               Teuchos::RCP<TreeVector> u) override;

  virtual AmanziSolvers::FnBaseDefs::ModifyCorrectionResult ModifyCorrection(
    double h,
    Teuchos::RCP<const TreeVector> res,
    Teuchos::RCP<const TreeVector> u,
    Teuchos::RCP<TreeVector> du) override;

  // computes the non-linear functional g = g(t,u,udot)
  virtual void FunctionalResidual(double t_old,
                                  double t_new,
                                  Teuchos::RCP<const TreeVector> u_old,
                                  Teuchos::RCP<TreeVector> u_new,
                                  Teuchos::RCP<TreeVector> g) override;

  // -- Commit any secondary (dependent) variables.
  virtual void CommitStep(double t_old, double t_new, const Tag& tag) override;
  virtual void FailStep(double t_old, double t_new, const Tag& tag) override;

 protected:
  Key snow_dens_key_;
  Key snow_age_key_;
  Key new_snow_key_;
  Key snow_source_key_;
  Key snow_death_rate_key_;
  Key area_frac_key_;

  double density_snow_max_;


 private:
  // factory registration
  static RegisteredPKFactory<ImplicitSubgrid> reg_;
};

} // namespace SurfaceBalance
} // namespace Amanzi

#endif
