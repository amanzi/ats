/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! An implicit PK for surface balance snow SWE conservation.


/*!

This is a balance PK whose conserved quantity is snow SWE.  The energy balance
comes in as it provides the energy needed to melt snow.  So source terms
include snow precipitation and snowmelt.  It also manages snow density, which
should get rethought a bit.

There is also some wierd hackiness here about area fractions -- see ATS Issue
#8

.. _subgrid-balance-pk-spec:
.. admonition:: subgrid-balance-pk-spec

    * `"absolute error tolerance`" ``[double]`` **0.01** ``[m]``    

    INCLUDES:

    - ``[balance-pk-spec]`` This *is a* `Balance Equation`_

    Not typically set by user, defaults work:

    * `"conserved quantity key`" ``[string]`` **LAYER-snow_water_equivalent**
      Sets the default conserved quantity key, so this is likely not supplied
      by the user. `[m]`
    * `"snow density key`" ``[string]`` **LAYER-density** Default snow density
      key. `[kg m^-3]`
    * `"snow age key`" ``[string]`` **LAYER-age** Default snow age key. `[d]`
    * `"new snow key`" ``[string]`` **LAYER-source** Default new snow key. `[m SWE s^-1]`
    * `"area fractions key`" ``[string]`` **LAYER-fractional_areas** Subgrid
      model fractional areas, see note above. `[-]`
    * `"snow death rate key`" ``[string]`` **LAYER-death_rate** Deals with last
      tiny bit of snowmelt.
    
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

  // main methods
  // -- Setup data.
  virtual void Setup(const Teuchos::Ptr<State>& S) override;

  // -- Initialize owned (dependent) variables.
  virtual void Initialize(const Teuchos::Ptr<State>& S) override;

  virtual bool ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0,
          Teuchos::RCP<TreeVector> u) override;


  virtual AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
  ModifyCorrection(double h, Teuchos::RCP<const TreeVector> res,
                   Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> du) override;
  
  // computes the non-linear functional g = g(t,u,udot)
  virtual void FunctionalResidual(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                   Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g) override;

  // -- Commit any secondary (dependent) variables.
  virtual void CommitStep(double t_old, double t_new,  const Teuchos::RCP<State>& S) override;

 protected:
  // multiple primary variables
  Teuchos::RCP<PrimaryVariableFieldEvaluator> pvfe_snow_dens_;
  Teuchos::RCP<PrimaryVariableFieldEvaluator> pvfe_snow_death_rate_;

  Key snow_dens_key_;
  Key snow_age_key_;
  Key new_snow_key_;
  Key snow_source_key_;
  Key snow_death_rate_key_;
  Key area_frac_key_;

 private:
  // factory registration
  static RegisteredPKFactory<ImplicitSubgrid> reg_;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
