/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------
   ATS

   DOCUMENT ME
   Surface Energy Balance for Snow Surface and Ground Surface
   Calculates Energy flux, rate or water, and water temperature
   entering through the surface skin.  Snow surface energy balance
   is calculated at equilibrium with ground/surface water and Air.

   ------------------------------------------------------------------------- */

#include <algorithm>

#include "pk_helpers.hh"
#include "seb_physics_defs.hh"
#include "surface_balance_implicit_subgrid.hh"

namespace Amanzi {
namespace SurfaceBalance {

ImplicitSubgrid::ImplicitSubgrid(Teuchos::ParameterList& pk_tree,
                                 const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                                 const Teuchos::RCP<State>& S,
                                 const Teuchos::RCP<TreeVector>& solution)
  : PK(pk_tree, global_list, S, solution), SurfaceBalanceBase(pk_tree, global_list, S, solution)
{
  if (!plist_->isParameter("conserved quantity key suffix"))
    plist_->set("conserved quantity key suffix", "snow_water_equivalent");

  // set up keys
  Key domain_surf = Keys::readDomainHint(*plist_, domain_, "snow", "surface");
  snow_dens_key_ = Keys::readKey(*plist_, domain_, "snow density", "density");
  snow_age_key_ = Keys::readKey(*plist_, domain_, "snow age", "age");
  new_snow_key_ = Keys::readKey(*plist_, domain_, "new snow source", "source");
  snow_death_rate_key_ = Keys::readKey(*plist_, domain_, "snow death rate", "death_rate");

  density_snow_max_ = plist_->get<double>("max density of snow [kg m^-3]", 600.);

  // set the error tolerance for snow
  plist_->set("absolute error tolerance", 0.01);
}

// main methods
// -- Setup data.
void
ImplicitSubgrid::Setup()
{
  SurfaceBalanceBase::Setup();

  // requirements: things I use
  requireAtNext(new_snow_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  // requirements: other primary variables
  requireAtNext(snow_dens_key_, tag_next_, *S_, name_)
    .SetMesh(mesh_)
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  requireAtCurrent(snow_dens_key_, tag_current_, *S_, name_);

  requireAtNext(snow_death_rate_key_, tag_next_, *S_, name_)
    .SetMesh(mesh_)
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  requireAtCurrent(snow_death_rate_key_, tag_current_, *S_, name_);

  requireAtNext(snow_age_key_, tag_next_, *S_, name_)
    .SetMesh(mesh_)
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  requireAtCurrent(snow_age_key_, tag_current_, *S_, name_);
}

// -- Initialize owned (dependent) variables.
void
ImplicitSubgrid::Initialize()
{
  SurfaceBalanceBase::Initialize();

  // initialize snow density, age
  Teuchos::ParameterList& ic_list = plist_->sublist("initial condition");

  if (!S_->GetRecord(snow_dens_key_, tag_next_).initialized()) {
    if (ic_list.isParameter("restart file")) {
      // initialize density, age from restart file
      S_->GetRecordW(snow_dens_key_, tag_next_, name_).Initialize(ic_list);
    } else if (plist_->isSublist("initial condition snow density")) {
      S_->GetRecordW(snow_dens_key_, tag_next_, name_)
        .Initialize(plist_->sublist("initial condition snow density"));
    } else {
      // initialize density to fresh powder, age to 0
      Relations::ModelParams params;
      S_->GetW<CompositeVector>(snow_dens_key_, tag_next_, name_)
        .PutScalar(params.density_freshsnow);
      S_->GetRecordW(snow_dens_key_, tag_next_, name_).set_initialized();
    }
  }

  if (!S_->GetRecord(snow_age_key_, tag_next_).initialized()) {
    if (ic_list.isParameter("restart file")) {
      // initialize density, age from restart file
      S_->GetRecordW(snow_age_key_, tag_next_, name_).Initialize(ic_list);
    } else if (plist_->isSublist("initial condition snow age")) {
      S_->GetRecordW(snow_age_key_, tag_next_, name_)
        .Initialize(plist_->sublist("initial condition snow age"));
    } else {
      // initialize age to fresh powder, age to 0
      S_->GetW<CompositeVector>(snow_age_key_, tag_next_, name_).PutScalar(0.);
      S_->GetRecordW(snow_age_key_, tag_next_, name_).set_initialized();
    }
  }

  S_->GetW<CompositeVector>(snow_death_rate_key_, tag_next_, name_).PutScalar(0.);
  S_->GetRecordW(snow_death_rate_key_, tag_next_, name_).set_initialized();
}


bool
ImplicitSubgrid::ModifyPredictor(double h,
                                 Teuchos::RCP<const TreeVector> u0,
                                 Teuchos::RCP<TreeVector> u)
{
  Epetra_MultiVector& u_vec = *u->Data()->ViewComponent("cell", false);
  for (int c = 0; c != u_vec.MyLength(); ++c) { u_vec[0][c] = std::max(0., u_vec[0][c]); }
  return true;
}


// -- Modify the correction.
AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
ImplicitSubgrid::ModifyCorrection(double h,
                                  Teuchos::RCP<const TreeVector> res,
                                  Teuchos::RCP<const TreeVector> u,
                                  Teuchos::RCP<TreeVector> du)
{
  Teuchos::OSTab tab = vo_->getOSTab();

  // modify correction to enforce nonnegativity
  int n_modified = 0;
  const Epetra_MultiVector& snow_depth = *u->Data()->ViewComponent("cell", false);
  Epetra_MultiVector& dsnow_depth = *du->Data()->ViewComponent("cell", false);
  for (int c = 0; c != snow_depth.MyLength(); ++c) {
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
ImplicitSubgrid::FunctionalResidual(double t_old,
                                    double t_new,
                                    Teuchos::RCP<TreeVector> u_old,
                                    Teuchos::RCP<TreeVector> u_new,
                                    Teuchos::RCP<TreeVector> g)
{
  int cycle = S_->get_cycle(tag_next_);

  // first calculate the "snow death rate", or rate of snow SWE that must melt over this
  // timestep if the snow is to go to zero.
  auto& snow_death_rate =
    *S_->GetW<CompositeVector>(snow_death_rate_key_, tag_next_, name_).ViewComponent("cell", false);
  S_->GetEvaluator(cell_vol_key_, tag_next_).Update(*S_, name_);
  const auto& cell_volume =
    *S_->Get<CompositeVector>(cell_vol_key_, tag_next_).ViewComponent("cell", false);
  snow_death_rate.PutScalar(0.);

  //S_->GetEvaluator(conserved_key_, tag_current_).Update(*S_, name_);
  S_->GetEvaluator(conserved_key_, tag_next_).Update(*S_, name_);
  const auto& swe_old_v =
    *S_->Get<CompositeVector>(conserved_key_, tag_current_).ViewComponent("cell", false);
  const auto& swe_new_v =
    *S_->Get<CompositeVector>(conserved_key_, tag_next_).ViewComponent("cell", false);
  for (int c = 0; c != snow_death_rate.MyLength(); ++c) {
    if (swe_new_v[0][c] <= 0.) {
      snow_death_rate[0][c] = swe_old_v[0][c] / (t_new - t_old) / cell_volume[0][c];
    }
  }
  // This line of code is commented out to ensure consistent code with master.
  // Uncommenting removes the bug described in ats#123, but does not fix it,
  // because it results in the cycle.
  // changedEvaluatorPrimary(snow_death_rate_key_, tag_next_, *S_);

  // update the residual
  SurfaceBalanceBase::FunctionalResidual(t_old, t_new, u_old, u_new, g);

  // now fill the role of age/density evaluator, as these depend upon old and new values
  const auto& snow_age_old =
    *S_->Get<CompositeVector>(snow_age_key_, tag_current_).ViewComponent("cell", false);
  auto& snow_age_new =
    *S_->GetW<CompositeVector>(snow_age_key_, tag_next_, name_).ViewComponent("cell", false);

  const auto& snow_dens_old =
    *S_->Get<CompositeVector>(snow_dens_key_, tag_current_).ViewComponent("cell", false);
  auto& snow_dens_new =
    *S_->GetW<CompositeVector>(snow_dens_key_, tag_next_, name_).ViewComponent("cell", false);

  S_->GetEvaluator(new_snow_key_, tag_next_).Update(*S_, name_);
  const auto& new_snow =
    *S_->Get<CompositeVector>(new_snow_key_, tag_next_).ViewComponent("cell", false);

  S_->GetEvaluator(source_key_, tag_next_).Update(*S_, name_);
  const auto& source =
    *S_->Get<CompositeVector>(source_key_, tag_next_).ViewComponent("cell", false);

  Relations::ModelParams params;
  double dt_days = (t_new - t_old) / 86400.;
  for (int c = 0; c != snow_dens_new.MyLength(); ++c) {
    double swe_added = new_snow[0][c] * (t_new - t_old) * cell_volume[0][c];
    double swe_lost = (new_snow[0][c] - source[0][c]) * (t_new - t_old) * cell_volume[0][c];
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
      snow_age_new[0][c] =
        (age_settled * std::max(swe_old - swe_lost, 0.) + age_new_snow * swe_added) /
        (std::max(swe_old - swe_lost, 0.) + swe_added);
      snow_dens_new[0][c] =
        (dens_settled * std::max(swe_old - swe_lost, 0.) + params.density_freshsnow * swe_added) /
        (std::max(swe_old - swe_lost, 0.) + swe_added);
      snow_dens_new[0][c] = std::min(snow_dens_new[0][c], density_snow_max_);
    }
  }

  // This line of code is commented out to ensure consistent code with master.
  // Uncommenting removes the bug described in ats#123, but does not fix it,
  // because it results in the cycle.
  // changedEvaluatorPrimary(snow_dens_key_, tag_next_, *S_);
  // changedEvaluatorPrimary(snow_age_key_, tag_next_, *S_);

  // debugging
  std::vector<std::string> vnames;
  vnames.push_back("snow age");
  vnames.push_back("snow dens");

  std::vector<Teuchos::Ptr<const CompositeVector>> vecs;
  vecs.push_back(S_->GetPtr<CompositeVector>(snow_age_key_, tag_next_).ptr());
  vecs.push_back(S_->GetPtr<CompositeVector>(snow_dens_key_, tag_next_).ptr());
  db_->WriteVectors(vnames, vecs, false);
}


void
ImplicitSubgrid::CommitStep(double t_old, double t_new, const Tag& tag_next)
{
  SurfaceBalanceBase::CommitStep(t_old, t_new, tag_next);

  AMANZI_ASSERT(tag_next == tag_next_ || tag_next == Tags::NEXT);
  Tag tag_current = tag_next == tag_next_ ? tag_current_ : Tags::CURRENT;

  assign(snow_age_key_, tag_current, tag_next, *S_);
  assign(snow_dens_key_, tag_current, tag_next, *S_);
  assign(snow_death_rate_key_, tag_current, tag_next, *S_);
}


void
ImplicitSubgrid::FailStep(double t_old, double t_new, const Tag& tag)
{
  SurfaceBalanceBase::FailStep(t_old, t_new, tag);
  if (tag == tag_next_) {
    S_->Assign(snow_age_key_, tag_next_, tag_current_);
    changedEvaluatorPrimary(snow_age_key_, tag_next_, *S_);

    S_->Assign(snow_dens_key_, tag_next_, tag_current_);
    changedEvaluatorPrimary(snow_dens_key_, tag_next_, *S_);

    S_->Assign(snow_death_rate_key_, tag_next_, tag_current_);
    changedEvaluatorPrimary(snow_death_rate_key_, tag_next_, *S_);
  }
}

} // namespace SurfaceBalance
} // namespace Amanzi
