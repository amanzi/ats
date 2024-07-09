/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Distributes and downregulates potential transpiration to the rooting zone.
#include "Brent.hh"
#include "transpiration_distribution_relperm_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

SoilPlantFluxFunctor::SoilPlantFluxFunctor(AmanziMesh::Entity_ID sc_,
                                           const AmanziMesh::MeshCache::cEntity_ID_View& cells_of_col_,
                                           const LandCover& lc_,
                                           const cView_type& soil_pc_,
                                           const cView_type& soil_kr_,
                                           const cView_type& f_root_,
                                           const cView_type& pet_,
                                           const cView_type& cv_,
                                           const cView_type& sa_,
                                           double c0_,
                                           double krp_,
                                           double rho_,
                                           double g_)
  : sc(sc_),
    cells_of_col(cells_of_col_),
    lc(lc_),
    soil_kr(soil_kr_),
    soil_pc(soil_pc_),
    f_root(f_root_),
    pet(pet_),
    cv(cv_),
    sa(sa_),
    c0(c0_),
    krp(krp_),
    rho_g(rho_ * g_)
{}


// Functor used in the solve for plant water potential
//
// Note the sign here -- pet should not change, as plant water potential
// increases, soil->plant flux should increase, so total_trans should increase,
// so this function is hopefully monotonically decreasing with plant_pc
double
SoilPlantFluxFunctor::operator()(double plant_pc) const
{
  double beta = computeTranspirationReductionFunction(plant_pc);
  double total_trans = 0.;

  // root pc is given by plant pc (assumed at the surface) - hydrostatic
  double root_pc = plant_pc;
  double sa_col = sa(sc, 0) * 2.0;
  for (auto c : cells_of_col) {
    double dz_on_2 = cv(c, 0) / sa_col;

    // top half-cell, note g is negative
    root_pc += rho_g * dz_on_2;

    // compute flux
    total_trans += computeSoilPlantFlux(root_pc, c);

    // bottom half-cell
    root_pc += rho_g * dz_on_2;
  }
  return beta * pet(sc, 0) - total_trans;
}


// beta in left hand side
double
SoilPlantFluxFunctor::computeTranspirationReductionFunction(double plant_pc) const
{
  if (plant_pc <= lc.stomata_open_capillary_pressure) {
    return 1.0;
  } else if (plant_pc >= lc.stomata_closed_capillary_pressure) {
    return 0.0;
  } else {
    return (plant_pc - lc.stomata_closed_capillary_pressure) /
           (lc.stomata_open_capillary_pressure - lc.stomata_closed_capillary_pressure);
  }
}


// right hand side
double
SoilPlantFluxFunctor::computeSoilPlantFlux(double root_pc, AmanziMesh::Entity_ID c) const
{
  double kr = (root_pc > soil_pc(c, 0)) ? soil_kr(c, 0) : krp;
  return -c0 * f_root(c, 0) * kr * (soil_pc(c, 0) - root_pc);
}


// all right hand sides
void
SoilPlantFluxFunctor::computeSoilPlantFluxes(double plant_pc, const View_type& res) const
{
  double root_pc = plant_pc;
  double sa_col = sa(sc, 0) * 2.0;
  for (auto c : cells_of_col) {
    double dz_on_2 = cv(c, 0) / sa_col;

    // top half-cell
    root_pc += rho_g * dz_on_2;

    // compute flux
    res(c, 0) = computeSoilPlantFlux(root_pc, c) * sa(sc, 0) / cv(c, 0);

    // bottom half-cell
    root_pc += rho_g * dz_on_2;
  }
}


const std::string TranspirationDistributionRelPermEvaluator::eval_type =
  "transpiration distribution, relative permeability";

// Constructor from ParameterList
TranspirationDistributionRelPermEvaluator::TranspirationDistributionRelPermEvaluator(
  const Teuchos::RCP<Teuchos::ParameterList>& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  InitializeFromPlist_();
}


// Virtual copy constructor
Teuchos::RCP<Evaluator>
TranspirationDistributionRelPermEvaluator::Clone() const
{
  return Teuchos::rcp(new TranspirationDistributionRelPermEvaluator(*this));
}


// Initialize by setting up dependencies
void
TranspirationDistributionRelPermEvaluator::InitializeFromPlist_()
{
  // Set up my dependencies
  // - defaults to prefixed via domain
  domain_sub_ = Keys::getDomain(my_keys_.front().first);
  domain_surf_ = Keys::readDomainHint(*plist_, domain_sub_, "domain", "surface");
  Tag tag = my_keys_.front().second;

  // tolerances for root finding
  tol_ = plist_->get<double>("tolerance", 1.e-12);
  nits_ = plist_->get<int>("maximum number of iterations", 100);

  // note this should probably not include salinity, etc?
  rho_ = plist_->get<double>("water density in plant [kg m^-3]", 1000.);

  // my_keys_.front() will always (for now?) be transpiration...
  Key plant_pc_key =
    Keys::readKey(*plist_, domain_surf_, "plant capillary pressure", "capillary_pressure_plant");
  my_keys_.emplace_back(std::make_pair(plant_pc_key, tag));

  // - pull Keys from plist
  // dependency: soil p
  soil_pc_key_ =
    Keys::readKey(*plist_, domain_sub_, "soil capillary pressure", "capillary_pressure_gas_liq");
  dependencies_.insert(KeyTag{ soil_pc_key_, tag });

  // dependency: soil rel perm
  soil_kr_key_ = Keys::readKey(
    *plist_, domain_sub_, "soil hydraulic conductivity", "relative_hydraulic_conductivity");
  dependencies_.insert(KeyTag{ soil_kr_key_, tag });

  // dependency: rooting_depth_fraction
  f_root_key_ = Keys::readKey(*plist_, domain_sub_, "plant rooting fraction", "root_fraction");
  dependencies_.insert(KeyTag{ f_root_key_, tag });

  // dependency: cv, sa
  cv_key_ = Keys::readKey(*plist_, domain_sub_, "cell volume", "cell_volume");
  dependencies_.insert(KeyTag{ cv_key_, tag });
  sa_key_ = Keys::readKey(*plist_, domain_surf_, "surface area", "cell_volume");
  dependencies_.insert(KeyTag{ sa_key_, tag });

  // dependency: transpiration
  potential_trans_key_ =
    Keys::readKey(*plist_, domain_surf_, "potential transpiration", "potential_transpiration_mols");
  dependencies_.insert(KeyTag{ potential_trans_key_, tag });

  // c0, maximal trans rate
  c0_ = plist_->get<double>("total maximal conductance [mol m^-2 s^-1 MPa^-1]", 10.) *
        1.e-6; // per MPa --> per Pa
  krp_ = plist_->get<double>("plant relative conductance [-]", 0.);

  land_cover_ =
    getLandCoverMap(plist_->sublist("model parameters"),
                    { "stomata_closed_capillary_pressure", "stomata_open_capillary_pressure" });
}


void
TranspirationDistributionRelPermEvaluator::Evaluate_(const State& S,
                                                     const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;

  // on the subsurface
  auto soil_pc = S.Get<CompositeVector>(soil_pc_key_, tag).viewComponent("cell", false);
  auto soil_kr = S.Get<CompositeVector>(soil_kr_key_, tag).viewComponent("cell", false);
  auto f_root = S.Get<CompositeVector>(f_root_key_, tag).viewComponent("cell", false);
  auto cv = S.Get<CompositeVector>(cv_key_, tag).viewComponent("cell", false);
  auto sa = S.Get<CompositeVector>(sa_key_, tag).viewComponent("cell", false);

  const auto& gravity = S.Get<AmanziGeometry::Point>("gravity", Tags::DEFAULT);
  double g = gravity[gravity.dim() - 1];

  // on the surface
  auto potential_trans =
    S.Get<CompositeVector>(potential_trans_key_, tag).viewComponent("cell", false);

  auto trans_v = result[0]->viewComponent("cell", false);
  auto plant_pc_v = result[1]->viewComponent("cell", false);

  auto surf_mesh = S.GetMesh(domain_surf_);
  const AmanziMesh::MeshCache& subsurf_mesh = S.GetMesh(domain_sub_)->getCache();

  for (const auto& region_lc : land_cover_) {
    auto lc_ids = surf_mesh->getSetEntities<MemSpace_kind::DEVICE>(
      region_lc.first, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
    const LandCover& lc_pars = region_lc.second;
    double krp(krp_), c0(c0_), rho(rho_), tol(tol_);

    if (lc_ids.size() > 0) {
      int nits(nits_);

      Kokkos::parallel_for(
        "TranspirationDistributionRelPermEvaluator::Evaluate",
        lc_ids.size(),
        KOKKOS_LAMBDA(const int i) {
          AmanziMesh::Entity_ID sc = lc_ids(i);
          if (potential_trans(sc, 0) > 0. || krp > 0.) {
            SoilPlantFluxFunctor func(sc,
                    subsurf_mesh.columns.getCells<MemSpace_kind::DEVICE>(sc),
                    lc_pars,
                    soil_pc,
                    soil_kr,
                    f_root,
                    potential_trans,
                    cv,
                    sa,
                    c0,
                    krp,
                    rho,
                    g);

            // bracket the root -- linear to the left of 1, log to the right of 1
            Kokkos::pair<double, double> ab;
            if (func(0.) > 0.) {
              ab.first = 0.;
              ab.second = 1.e4;
              while (func(ab.second) > 0) {
                ab.first = ab.second;
                ab.second *= 10;
                assert(ab.second < 1.e14); // failed to bracket root
              }
            } else {
              ab.second = 0.;
              ab.first = -1.e4;
              while (func(ab.first) < 0) {
                ab.second = ab.first;
                ab.first *= 10;
                assert(ab.first > -1.e14); // failed to bracket root
              }
            }

            // compute the plant capillary pressure using a root-finder
            int itrs = nits;
            plant_pc_v(sc, 0) =
              Amanzi::Utils::findRootBrent(func, ab.first, ab.second, tol, &itrs);
            assert(itrs > 0 && itrs <= nits);

            // compute the distributed transpiration fluxes for each grid cell
            func.computeSoilPlantFluxes(plant_pc_v(sc, 0), trans_v);

          } else {
            for (auto c : subsurf_mesh.columns.getCells<MemSpace_kind::DEVICE>(sc)) {
              trans_v(c, 0) = 0.;
            }
          }
        });
    }
  }
  AMANZI_ASSERT(!std::isnan(trans_v(0, 0)));
}


void
TranspirationDistributionRelPermEvaluator::EvaluatePartialDerivative_(
  const State& S,
  const Key& wrt_key,
  const Tag& wrt_tag,
  const std::vector<CompositeVector*>& result)
{
  AMANZI_ASSERT(false);
}


void
TranspirationDistributionRelPermEvaluator::EnsureCompatibility_ToDeps_(State& S)
{
  Tag tag = my_keys_.front().second;

  // Create an unowned factory to check my dependencies.
  // -- first those on the subsurface mesh
  CompositeVectorSpace dep_fac;
  dep_fac.SetMesh(S.GetMesh(domain_sub_))->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  S.Require<CompositeVector, CompositeVectorSpace>(f_root_key_, tag).Update(dep_fac);
  S.Require<CompositeVector, CompositeVectorSpace>(soil_pc_key_, tag).Update(dep_fac);
  S.Require<CompositeVector, CompositeVectorSpace>(soil_kr_key_, tag).Update(dep_fac);
  S.Require<CompositeVector, CompositeVectorSpace>(cv_key_, tag).Update(dep_fac);

  // -- next those on the surface mesh
  CompositeVectorSpace surf_fac;
  surf_fac.SetMesh(S.GetMesh(domain_surf_))->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  S.Require<CompositeVector, CompositeVectorSpace>(potential_trans_key_, tag).Update(surf_fac);
  S.Require<CompositeVector, CompositeVectorSpace>(sa_key_, tag).Update(surf_fac);

  // gravity
  S.Require<AmanziGeometry::Point>("gravity", Tags::NEXT);
}


void
TranspirationDistributionRelPermEvaluator::EnsureCompatibility_Structure_(State& S)
{
  Tag tag = my_keys_.front().second;

  // my_keys_[0] is transpiration
  S.Require<CompositeVector, CompositeVectorSpace>(my_keys_.front().first, tag)
    .SetMesh(S.GetMesh(domain_sub_))
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  // my_keys_[1] is plant water potential
  S.Require<CompositeVector, CompositeVectorSpace>(my_keys_.back().first, tag)
    .SetMesh(S.GetMesh(domain_surf_))
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
}


} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
