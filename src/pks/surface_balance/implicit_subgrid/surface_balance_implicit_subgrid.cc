/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon

   DOCUMENT ME
   Surface Energy Balance for Snow Surface and Ground Surface
   Calculates Energy flux, rate or water, and water temperature
   entering through the surface skin.  Snow surface energy balance
   is calculated at equilibrium with ground/surface water and Air.

   ------------------------------------------------------------------------- */

#include <algorithm>
#include "boost/algorithm/string/predicate.hpp"

#include "seb_physics_defs.hh"
#include "seb_physics_funcs.hh"
#include "surface_balance_implicit_subgrid.hh"

namespace Amanzi {
namespace SurfaceBalance {

ImplicitSubgrid::ImplicitSubgrid(Teuchos::ParameterList& pk_tree,
                                     const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                                     const Teuchos::RCP<State>& S,
                                     const Teuchos::RCP<TreeVector>& solution):
  PK(pk_tree, global_list,  S, solution),
  SurfaceBalanceBase(pk_tree, global_list,  S, solution)
{
  if (!plist_->isParameter("conserved quantity key suffix"))
    plist_->set("conserved quantity key suffix", "snow_water_equivalent");

  Teuchos::ParameterList& FElist = S->FEList();

  // set up keys
  Key domain_surf;
  if (domain_ == "snow") {
    domain_surf = "surface";
  } else if (boost::starts_with(domain_, "snow_")) {
    domain_surf = std::string("surface_") + domain_.substr(5,domain_.size());
  } else {
    Errors::Message message("SurfaceBalance::ImplicitSubgrid PK cannot deduce surface domain name.");
    Exceptions::amanzi_throw(message);
  }
  
  snow_dens_key_ = Keys::readKey(*plist_, domain_, "snow density", "density");
  snow_age_key_ = Keys::readKey(*plist_, domain_, "snow age", "age");
  new_snow_key_ = Keys::readKey(*plist_, domain_, "new snow source", "source");
  area_frac_key_ = Keys::readKey(*plist_, domain_surf, "area fractions", "fractional_areas");
  snow_death_rate_key_ = Keys::readKey(*plist_, domain_, "snow death rate", "death_rate");
  
  // set up additional primary variables -- this is very hacky, and can become an evaluator in new-state
  // -- snow density
  Teuchos::ParameterList& snow_dens_sublist = FElist.sublist(snow_dens_key_);
  snow_dens_sublist.set("evaluator name", snow_dens_key_);
  snow_dens_sublist.set("field evaluator type", "primary variable");
  S->FEList().set(snow_dens_key_, snow_dens_sublist);

  // -- snow death rate
  Teuchos::ParameterList& snow_death_rate_sublist = FElist.sublist(snow_death_rate_key_);
  snow_death_rate_sublist.set("evaluator name", snow_death_rate_key_);
  snow_death_rate_sublist.set("field evaluator type", "primary variable");
  S->FEList().set(snow_death_rate_key_, snow_death_rate_sublist);
  
  // set the error tolerance for snow
  plist_->set("absolute error tolerance", 0.01);
}

// main methods
// -- Setup data.
void
ImplicitSubgrid::Setup(const Teuchos::Ptr<State>& S) {
  SurfaceBalanceBase::Setup(S);

  // requireiments: things I use
  S->RequireField(new_snow_key_)->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator(new_snow_key_);
  
  // requirements: other primary variables
  Teuchos::RCP<FieldEvaluator> fm;
  S->RequireField(snow_dens_key_, name_)->SetMesh(mesh_)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator(snow_dens_key_);
  fm = S->GetFieldEvaluator(snow_dens_key_);
  pvfe_snow_dens_ = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(fm);
  if (pvfe_snow_dens_ == Teuchos::null) {
    Errors::Message message("SurfaceBalanceSEB: error, failure to initialize primary variable");
    Exceptions::amanzi_throw(message);
  }

  S->RequireField(snow_death_rate_key_, name_)->SetMesh(mesh_)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator(snow_death_rate_key_);
  fm = S->GetFieldEvaluator(snow_death_rate_key_);
  pvfe_snow_death_rate_ = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(fm);
  if (pvfe_snow_death_rate_ == Teuchos::null) {
    Errors::Message message("SurfaceBalanceSEB: error, failure to initialize primary variable");
    Exceptions::amanzi_throw(message);
  }
  
  // requirements: internally we must track snow age, this should become an
  // evaluator with snow_density.
  S->RequireField(snow_age_key_, name_)->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  // requirements: area fractions will be updated by us to make sure the old
  // time's value is used.  It is likely an error if anything else depends
  // upon area fractions except for the seb_evaluator (which explicitly
  // comments out the dependency) and us (which manage when we update it to
  // not do it until commit state when we update to the new value.
  // :ISSUE:#8
  S->RequireFieldEvaluator(area_frac_key_);
}

// -- Initialize owned (dependent) variables.
void
ImplicitSubgrid::Initialize(const Teuchos::Ptr<State>& S) {
  SurfaceBalanceBase::Initialize(S);

  // initialize snow density, age
  AMANZI_ASSERT(plist_->isSublist("initial condition"));
  Teuchos::ParameterList& ic_list = plist_->sublist("initial condition");
  
  if (!S->GetField(snow_dens_key_)->initialized()) {
    if (ic_list.isParameter("restart file")) {
      // initialize density, age from restart file
      S->GetField(snow_dens_key_, name_)->Initialize(ic_list);
      S->GetField(snow_dens_key_, name_)->set_initialized();
    } else if (plist_->isSublist("initial condition snow density")) {
      S->GetField(snow_dens_key_, name_)->Initialize(plist_->sublist("initial condition snow density"));
    } else {
      // initialize density to fresh powder, age to 0
      SEBPhysics::ModelParams params;
      S->GetFieldData(snow_dens_key_,name_)->PutScalar(params.density_freshsnow);
      S->GetField(snow_dens_key_, name_)->set_initialized();
    }
  }

  if (!S->GetField(snow_age_key_)->initialized()) {
    if (ic_list.isParameter("restart file")) {
      // initialize density, age from restart file
      S->GetField(snow_age_key_, name_)->Initialize(ic_list);
      S->GetField(snow_age_key_, name_)->set_initialized();
    } else if (plist_->isSublist("initial condition snow age")) {
      S->GetField(snow_age_key_, name_)->Initialize(plist_->sublist("initial condition snow age"));
    } else {
      // initialize age to fresh powder, age to 0
      S->GetFieldData(snow_age_key_,name_)->PutScalar(0.);
      S->GetField(snow_age_key_, name_)->set_initialized();
    }
  }

  S->GetFieldData(snow_death_rate_key_,name_)->PutScalar(0.);
  S->GetField(snow_death_rate_key_, name_)->set_initialized();
}


bool
ImplicitSubgrid::ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0,
        Teuchos::RCP<TreeVector> u) {
  Epetra_MultiVector& u_vec = *u->Data()->ViewComponent("cell",false);
  for (int c=0; c!=u_vec.MyLength(); ++c) {
    u_vec[0][c] = std::max(0., u_vec[0][c]);
  }
  return true;
}


// -- Modify the correction.
AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
ImplicitSubgrid::ModifyCorrection(double h, Teuchos::RCP<const TreeVector> res,
        Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> du) {
  Teuchos::OSTab tab = vo_->getOSTab();

  // modify correction to enforce nonnegativity
  int n_modified = 0;
  const Epetra_MultiVector& snow_depth = *u->Data()->ViewComponent("cell",false);
  Epetra_MultiVector& dsnow_depth = *du->Data()->ViewComponent("cell",false);
  for (int c=0; c!=snow_depth.MyLength(); ++c) {
    if (snow_depth[0][c] - dsnow_depth[0][c] < 0.) {
      dsnow_depth[0][c] = snow_depth[0][c];
      n_modified++;
    }
  }

  // -- accumulate globally
  //  int n_modified_l = n_modified;
  //  u->Data()->Comm().SumAll(&n_modified_l, &n_modified, 1);
  return AmanziSolvers::FnBaseDefs::CORRECTION_NOT_MODIFIED; // ok to backtrack on this
}



// computes the non-linear functional g = g(t,u,udot)
void
ImplicitSubgrid::FunctionalResidual(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
        Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g) {
  int cycle = S_next_->cycle();

  AMANZI_ASSERT(S_next_->GetFieldEvaluator(new_snow_key_).get() == S_next_->GetFieldEvaluator(source_key_).get());

  // first calculate the "snow death rate", or rate of snow SWE that must melt over this
  // timestep if the snow is to go to zero.
  auto& snow_death_rate = *S_next_->GetFieldData(snow_death_rate_key_, name_)->ViewComponent("cell",false);
  S_next_->GetFieldEvaluator(cell_vol_key_)->HasFieldChanged(S_next_.ptr(), name_);
  const auto& cell_volume = *S_next_->GetFieldData(cell_vol_key_)->ViewComponent("cell",false);
  snow_death_rate.PutScalar(0.);

  S_inter_->GetFieldEvaluator(conserved_key_)->HasFieldChanged(S_inter_.ptr(), name_);
  S_next_->GetFieldEvaluator(conserved_key_)->HasFieldChanged(S_next_.ptr(), name_);
  const auto& swe_old_v = *S_inter_->GetFieldData(conserved_key_)->ViewComponent("cell", false);
  const auto& swe_new_v = *S_next_->GetFieldData(conserved_key_)->ViewComponent("cell", false);
  for (int c=0; c!=snow_death_rate.MyLength(); ++c) {
    if (swe_new_v[0][c] <= 0.) {
      snow_death_rate[0][c] = swe_old_v[0][c] / (t_new - t_old) / cell_volume[0][c];
    }
  }
  pvfe_snow_death_rate_->SetFieldAsChanged(S_next_.ptr());
  
  // update the residual
  SurfaceBalanceBase::FunctionalResidual(t_old, t_new, u_old, u_new, g);  

  // now fill the role of age/density evaluator, as these depend upon old and new values
  const auto& cell_vol = *S_next_->GetFieldData(cell_vol_key_)->ViewComponent("cell",false);

  const auto& snow_age_old = *S_inter_->GetFieldData(snow_age_key_)->ViewComponent("cell",false);
  auto& snow_age_new = *S_next_->GetFieldData(snow_age_key_, name_)->ViewComponent("cell",false);

  const auto& snow_dens_old = *S_inter_->GetFieldData(snow_dens_key_)->ViewComponent("cell",false);
  auto& snow_dens_new = *S_next_->GetFieldData(snow_dens_key_, name_)->ViewComponent("cell",false);
  
  S_next_->GetFieldEvaluator(new_snow_key_)->HasFieldChanged(S_next_.ptr(), name_);
  const auto& new_snow = *S_next_->GetFieldData(new_snow_key_)->ViewComponent("cell",false);

  S_next_->GetFieldEvaluator(source_key_)->HasFieldChanged(S_next_.ptr(), name_);
  const auto& source = *S_next_->GetFieldData(source_key_)->ViewComponent("cell",false);

  SEBPhysics::ModelParams params;
  double dt_days = (t_new - t_old) / 86400.;
  for (int c=0; c!=snow_dens_new.MyLength(); ++c) {
    double swe_added = new_snow[0][c] * (t_new - t_old) * cell_vol[0][c];
    double swe_lost = (new_snow[0][c] - source[0][c] ) * (t_new - t_old) * cell_vol[0][c];
    double swe_old = swe_old_v[0][c];

    double age_new_snow = dt_days / 2.;
    if (swe_old + swe_added - swe_lost < 1.e-10) {
      snow_age_new[0][c] = 0.;
      snow_dens_new[0][c] = params.density_freshsnow;
    } else {
      // age the old snow
      double age_settled = snow_age_old[0][c] + dt_days;
      double dens_settled = params.density_freshsnow * std::max(std::pow(age_settled, 0.3), 1.);
      
      // Match frost age with assigned density -- Calculate which day frost
      // density matched snow defermation function from (Martinec, 1977)
      //double age_frost = std::pow((params.density_frost / params.density_freshsnow),
      //        (1/0.3)) - 1 + dt;

      // ignoring frost, just weighting precip and old snow
      snow_age_new[0][c] = (age_settled * std::max(swe_old - swe_lost,0.) + age_new_snow * swe_added)
                           / (std::max(swe_old - swe_lost,0.) + swe_added);
      snow_dens_new[0][c] = (dens_settled * std::max(swe_old - swe_lost,0.) + params.density_freshsnow * swe_added)
                            / (std::max(swe_old - swe_lost,0.) + swe_added);
      snow_dens_new[0][c] = std::min(snow_dens_new[0][c], params.density_snow_max);
    }
  }
  pvfe_snow_dens_->SetFieldAsChanged(S_next_.ptr());

  // debugging
  std::vector<std::string> vnames;
  vnames.push_back("snow age");
  vnames.push_back("snow dens");

  std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
  vecs.push_back(S_next_->GetFieldData(snow_age_key_).ptr());
  vecs.push_back(S_next_->GetFieldData(snow_dens_key_).ptr());
  db_->WriteVectors(vnames, vecs, false);
}

void
ImplicitSubgrid::CommitStep(double t_old, double t_new,  const Teuchos::RCP<State>& S)
{
  // now update area frac
  S->GetFieldEvaluator(area_frac_key_)->HasFieldChanged(S.ptr(), name_);
  SurfaceBalanceBase::CommitStep(t_old, t_new, S);
}

} // namespace
} // namespace
