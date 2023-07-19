/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon, Chonggang Xu
*/

/* -------------------------------------------------------------------------
   ATS

   Simple implementation of CLM's Century model for carbon decomposition and a
   simplified 2-PFT (sedge, moss) vegetation model for creating carbon.

   CURRENT ASSUMPTIONS:
     1. parallel decomp not in the vertical
     2. fields are not ordered along the column, and so must be copied
     3. all columns have the same number of cells
   ------------------------------------------------------------------------- */

#include "MeshPartition.hh"
#include "pk_helpers.hh"
#include "bgc_simple_funcs.hh"

#include "bgc_simple.hh"

namespace Amanzi {
namespace BGC {


BGCSimple::BGCSimple(Teuchos::ParameterList& pk_tree,
                     const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                     const Teuchos::RCP<State>& S,
                     const Teuchos::RCP<TreeVector>& solution)
  : PK_Physical_Default(pk_tree, global_list, S, solution),
    PK(pk_tree, global_list, S, solution),
    ncells_per_col_(-1)
{
  // set up additional primary variables -- this is very hacky...
  // -- surface energy source
  domain_surf_ = Keys::readDomainHint(*plist_, domain_, "subsurface", "surface");

  // additional primary variables
  // -- transpiration
  trans_key_ = Keys::readKey(*plist_, domain_, "transpiration", "transpiration");
  // -- shortwave incoming shading
  shaded_sw_key_ = Keys::readKey(
    *plist_, domain_surf_, "shaded shortwave radiation", "shaded_shortwave_radiation");
  // -- lai
  total_lai_key_ =
    Keys::readKey(*plist_, domain_surf_, "total leaf area index", "total_leaf_area_index");

  // initial timestep
  dt_ = plist_->get<double>("initial time step", 1.);

  // Create the additional, non-managed data structures
  num_pools_ = plist_->get<int>("number of carbon pools", 7);

  // parameters
  lat_ = plist_->get<double>("latitude [degrees]");
  wind_speed_ref_ht_ = plist_->get<double>("wind speed reference height [m]", 2.0);
  cryoturbation_coef_ = plist_->get<double>("cryoturbation mixing coefficient [cm^2/yr]", 5.0);
  cryoturbation_coef_ /= 365.25e4; // convert to m^2/day
}

// is a PK
// -- Setup data
void
BGCSimple::Setup()
{
  PK_Physical_Default::Setup();

  // my mesh is the subsurface mesh, but we need the surface mesh, index by column, as well
  mesh_surf_ = S_->GetMesh(domain_surf_);

  // -- SoilCarbonParameters
  Teuchos::ParameterList& sc_params = plist_->sublist("soil carbon parameters");
  std::string mesh_part_name = sc_params.get<std::string>("mesh partition");
  const Functions::MeshPartition& mp = *S_->GetMeshPartition(mesh_part_name);
  const std::vector<std::string>& regions = mp.regions();

  for (const auto& region : regions) {
    sc_params_.push_back(
      Teuchos::rcp(new SoilCarbonParameters(num_pools_, sc_params.sublist(region))));
  }

  // -- PFTs -- old and new!
  Teuchos::ParameterList& pft_params = plist_->sublist("pft parameters");
  std::vector<std::string> pft_names;
  for (Teuchos::ParameterList::ConstIterator lcv = pft_params.begin(); lcv != pft_params.end();
       ++lcv) {
    std::string pft_name = lcv->first;
    pft_names.push_back(pft_name);
  }

  // set sizes
  num_pfts_ = pft_names.size();
  num_cols_ = mesh_surf_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  pfts_old_.resize(num_cols_);
  pfts_.resize(num_cols_);
  for (unsigned int col = 0; col != num_cols_; ++col) {
    int f = mesh_surf_->getEntityParent(AmanziMesh::Entity_kind::CELL, col);
    auto col_iter = mesh_->columns.getCells(col);
    std::size_t ncol_cells = col_iter.size();

    // unclear which this should be:
    // -- col area is the true face area
    double col_area = mesh_->getFaceArea(f);
    // -- col area is the projected face area
    // double col_area = mesh_surf_->getCellVolume(col);

    if (ncells_per_col_ < 0) {
      ncells_per_col_ = ncol_cells;
    } else {
      AMANZI_ASSERT(ncol_cells == ncells_per_col_);
    }

    pfts_old_[col].resize(num_pfts_);
    pfts_[col].resize(num_pfts_);

    for (int i = 0; i != num_pfts_; ++i) {
      std::string pft_name = pft_names[i];
      Teuchos::ParameterList& pft_plist = pft_params.sublist(pft_name);
      pfts_old_[col][i] = Teuchos::rcp(new PFT(pft_name, ncol_cells));
      pfts_old_[col][i]->Init(pft_plist, col_area);
      pfts_[col][i] = Teuchos::rcp(new PFT(*pfts_old_[col][i]));
    }
  }

  // -- soil carbon pools
  soil_carbon_pools_.resize(num_cols_);
  for (unsigned int col = 0; col != num_cols_; ++col) {
    soil_carbon_pools_[col].resize(ncells_per_col_);

    auto col_iter = mesh_->columns.getCells(col);
    ncells_per_col_ = col_iter.size();

    for (std::size_t i = 0; i != col_iter.size(); ++i) {
      // col_iter[i] = cell id, mp[cell_id] = index into partition list, sc_params_[index] = correct params
      soil_carbon_pools_[col][i] = Teuchos::rcp(new SoilCarbon(sc_params_[mp[col_iter[i]]]));
    }
  }

  // requirements: primary variable
  S_->Require<CompositeVector, CompositeVectorSpace>(key_, tag_next_, name_)
    .SetMesh(mesh_)
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, num_pools_);

  // requirements: other primary variables
  S_->Require<CompositeVector, CompositeVectorSpace>(trans_key_, tag_next_, name_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  requireEvaluatorPrimary(trans_key_, tag_next_, *S_);

  S_->Require<CompositeVector, CompositeVectorSpace>(shaded_sw_key_, tag_next_, name_)
    .SetMesh(mesh_surf_)
    ->SetGhosted()
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  requireEvaluatorPrimary(shaded_sw_key_, tag_next_, *S_);

  S_->Require<CompositeVector, CompositeVectorSpace>(total_lai_key_, tag_next_, name_)
    .SetMesh(mesh_surf_)
    ->SetGhosted()
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  requireEvaluatorPrimary(total_lai_key_, tag_next_, *S_);

  // requirement: diagnostics
  S_->Require<CompositeVector, CompositeVectorSpace>("co2_decomposition", tag_next_, name_)
    .SetMesh(mesh_)
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  S_->Require<CompositeVector, CompositeVectorSpace>("surface-total_biomass", tag_next_, name_)
    .SetMesh(mesh_surf_)
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, num_pfts_);
  S_->GetRecordSetW("surface-total_biomass").set_subfieldnames(pft_names);

  S_->Require<CompositeVector, CompositeVectorSpace>("surface-leaf_biomass", tag_next_, name_)
    .SetMesh(mesh_surf_)
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, num_pfts_);
  S_->GetRecordSetW("surface-leaf_biomass").set_subfieldnames(pft_names);

  S_->Require<CompositeVector, CompositeVectorSpace>("surface-leaf_area_index", tag_next_, name_)
    .SetMesh(mesh_surf_)
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, num_pfts_);
  S_->GetRecordSetW("surface-leaf_area_index").set_subfieldnames(pft_names);

  S_->Require<CompositeVector, CompositeVectorSpace>("surface-c_sink_limit", tag_next_, name_)
    .SetMesh(mesh_surf_)
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, num_pfts_);
  S_->GetRecordSetW("surface-c_sink_limit").set_subfieldnames(pft_names);

  S_->Require<CompositeVector, CompositeVectorSpace>(
      "surface-veg_total_transpiration", tag_next_, name_)
    .SetMesh(mesh_surf_)
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, num_pfts_);
  S_->GetRecordSetW("surface-veg_total_transpiration").set_subfieldnames(pft_names);

  // requirement: temp of each cell
  S_->RequireEvaluator("temperature", tag_next_);
  S_->Require<CompositeVector, CompositeVectorSpace>("temperature", tag_next_)
    .SetMesh(mesh_)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  // requirement: pressure
  S_->RequireEvaluator("pressure", tag_next_);
  S_->Require<CompositeVector, CompositeVectorSpace>("pressure", tag_next_)
    .SetMesh(mesh_)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  // requirements: surface cell volume
  S_->Require<CompositeVector, CompositeVectorSpace>("surface-cell_volume", tag_next_)
    .SetMesh(mesh_surf_)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  S_->RequireEvaluator("surface-cell_volume", tag_next_);

  // requirements: Met data
  S_->RequireEvaluator("surface-incoming_shortwave_radiation", tag_next_);
  S_->Require<CompositeVector, CompositeVectorSpace>("surface-incoming_shortwave_radiation",
                                                     tag_next_)
    .SetMesh(mesh_surf_)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  S_->RequireEvaluator("surface-air_temperature", tag_next_);
  S_->Require<CompositeVector, CompositeVectorSpace>("surface-air_temperature", tag_next_)
    .SetMesh(mesh_surf_)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  S_->RequireEvaluator("surface-vapor_pressure_air", tag_next_);
  S_->Require<CompositeVector, CompositeVectorSpace>("surface-vapor_pressure_air", tag_next_)
    .SetMesh(mesh_surf_)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  S_->RequireEvaluator("surface-wind_speed", tag_next_);
  S_->Require<CompositeVector, CompositeVectorSpace>("surface-wind_speed", tag_next_)
    .SetMesh(mesh_surf_)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  S_->RequireEvaluator("surface-co2_concentration", tag_next_);
  S_->Require<CompositeVector, CompositeVectorSpace>("surface-co2_concentration", tag_next_)
    .SetMesh(mesh_surf_)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
}

// -- Initialize owned (dependent) variables.
void
BGCSimple::Initialize()
{
  PK_Physical_Default::Initialize();

  // diagnostic variable
  S_->GetW<CompositeVector>("co2_decomposition", tag_next_, name_).PutScalar(0.);
  S_->GetRecordW("co2_decomposition", tag_next_, name_).set_initialized();

  S_->GetW<CompositeVector>("surface-c_sink_limit", tag_next_, name_).PutScalar(0.);
  S_->GetRecordW("surface-c_sink_limit", tag_next_, name_).set_initialized();

  S_->GetW<CompositeVector>(trans_key_, tag_next_, name_).PutScalar(0.);
  S_->GetRecordW(trans_key_, tag_next_, name_).set_initialized();

  S_->GetW<CompositeVector>("surface-total_biomass", tag_next_, name_).PutScalar(0.);
  S_->GetRecordW("surface-total_biomass", tag_next_, name_).set_initialized();

  S_->GetW<CompositeVector>("surface-leaf_area_index", tag_next_, name_).PutScalar(0.);
  S_->GetRecordW("surface-leaf_area_index", tag_next_, name_).set_initialized();

  S_->GetW<CompositeVector>("surface-veg_total_transpiration", tag_next_, name_).PutScalar(0.);
  S_->GetRecordW("surface-veg_total_transpiration", tag_next_, name_).set_initialized();

  // potentially initial aboveground vegetation data
  auto& leaf_biomass_field = S_->GetRecordW("surface-leaf_biomass", tag_next_, name_);
  if (!leaf_biomass_field.initialized()) {
    // -- Calculate the IC.
    if (plist_->isSublist("leaf biomass initial condition")) {
      Teuchos::ParameterList ic_plist = plist_->sublist("leaf biomass initial condition");
      leaf_biomass_field.Initialize(ic_plist);

      if (leaf_biomass_field.initialized()) {
        // -- copy into PFTs
        Epetra_MultiVector& bio =
          *S_->GetPtrW<CompositeVector>("surface-leaf_biomass", tag_next_, name_)
             ->ViewComponent("cell", false);

        int num_cols_ =
          mesh_surf_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
        for (int col = 0; col != num_cols_; ++col) {
          for (int i = 0; i != bio.NumVectors(); ++i) { pfts_old_[col][i]->Bleaf = bio[i][col]; }
        }
      }
    }

    if (!leaf_biomass_field.initialized()) {
      S_->GetW<CompositeVector>("surface-leaf_biomass", tag_next_, name_).PutScalar(0.);
      leaf_biomass_field.set_initialized();
    }
  }

  // init root carbon
  auto col_temp = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_depth = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_dz = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));

  S_->GetEvaluator("temperature", tag_next_).Update(*S_, name_);
  const Epetra_Vector& temp =
    *(*S_->Get<CompositeVector>("temperature", tag_next_).ViewComponent("cell", false))(0);

  int num_cols_ = mesh_surf_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  for (int col = 0; col != num_cols_; ++col) {
    FieldToColumn_(col, temp, col_temp.ptr());
    ColDepthDz_(col, col_depth.ptr(), col_dz.ptr());

    for (int i = 0; i != num_pfts_; ++i) {
      pfts_old_[col][i]->InitRoots(*col_temp, *col_depth, *col_dz);
    }
  }

  // ensure all initialization in both PFTs?  Not sure this is
  // necessary -- likely done in initial call to commit-state --etc
  for (int col = 0; col != num_cols_; ++col) {
    for (int i = 0; i != num_pfts_; ++i) { *pfts_[col][i] = *pfts_old_[col][i]; }
  }
}


// -- Commit any secondary (dependent) variables.
void
BGCSimple::CommitStep(double told, double tnew, const Tag& tag)
{
  // Copy the PFT over, which includes all additional state required, commit
  // the step as succesful.
  double dt = tnew - told;

  int num_cols_ = mesh_surf_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  for (int col = 0; col != num_cols_; ++col) {
    for (int i = 0; i != num_pfts_; ++i) { *pfts_old_[col][i] = *pfts_[col][i]; }
  }
}

// -- advance the model
bool
BGCSimple::AdvanceStep(double t_old, double t_new, bool reinit)
{
  double dt = t_new - t_old;

  Teuchos::OSTab out = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "----------------------------------------------------------------" << std::endl
               << "Advancing: t0 = " << S_->get_time(tag_current_)
               << " t1 = " << S_->get_time(tag_next_) << " h = " << dt << std::endl
               << "----------------------------------------------------------------" << std::endl;

  // Copy the PFT from old to new, in case we failed the previous attempt at
  // this timestep.  This is hackery to get around the fact that PFTs are not
  // (but should be) in state.
  AmanziMesh::Entity_ID num_cols_ =
    mesh_surf_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  for (AmanziMesh::Entity_ID col = 0; col != num_cols_; ++col) {
    for (int i = 0; i != num_pfts_; ++i) { *pfts_[col][i] = *pfts_old_[col][i]; }
  }

  // grab the required fields
  Epetra_MultiVector& sc_pools =
    *S_->GetW<CompositeVector>(key_, tag_next_, name_).ViewComponent("cell", false);
  Epetra_MultiVector& co2_decomp =
    *S_->GetW<CompositeVector>("co2_decomposition", tag_next_, name_).ViewComponent("cell", false);
  Epetra_MultiVector& trans =
    *S_->GetW<CompositeVector>(trans_key_, tag_next_, name_).ViewComponent("cell", false);
  Epetra_MultiVector& sw =
    *S_->GetW<CompositeVector>(shaded_sw_key_, tag_next_, name_).ViewComponent("cell", false);
  Epetra_MultiVector& biomass =
    *S_->GetW<CompositeVector>("surface-total_biomass", tag_next_, name_)
       .ViewComponent("cell", false);
  Epetra_MultiVector& leafbiomass =
    *S_->GetW<CompositeVector>("surface-leaf_biomass", tag_next_, name_)
       .ViewComponent("cell", false);
  Epetra_MultiVector& csink = *S_->GetW<CompositeVector>("surface-c_sink_limit", tag_next_, name_)
                                 .ViewComponent("cell", false);
  Epetra_MultiVector& total_lai =
    *S_->GetW<CompositeVector>(total_lai_key_, tag_next_, name_).ViewComponent("cell", false);
  Epetra_MultiVector& lai = *S_->GetW<CompositeVector>("surface-leaf_area_index", tag_next_, name_)
                               .ViewComponent("cell", false);
  Epetra_MultiVector& total_transpiration =
    *S_->GetW<CompositeVector>("surface-veg_total_transpiration", tag_next_, name_)
       .ViewComponent("cell", false);

  S_->GetEvaluator("temperature", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& temp =
    *S_->Get<CompositeVector>("temperature", tag_next_).ViewComponent("cell", false);

  S_->GetEvaluator("pressure", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& pres =
    *S_->Get<CompositeVector>("pressure", tag_next_).ViewComponent("cell", false);

  S_->GetEvaluator("surface-incoming_shortwave_radiation", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& qSWin =
    *S_->Get<CompositeVector>("surface-incoming_shortwave_radiation", tag_next_)
       .ViewComponent("cell", false);

  S_->GetEvaluator("surface-air_temperature", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& air_temp =
    *S_->Get<CompositeVector>("surface-air_temperature", tag_next_).ViewComponent("cell", false);

  S_->GetEvaluator("surface-vapor_pressure_air", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& vp_air =
    *S_->Get<CompositeVector>("surface-vapor_pressure_air", tag_next_).ViewComponent("cell", false);

  S_->GetEvaluator("surface-wind_speed", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& wind_speed =
    *S_->Get<CompositeVector>("surface-wind_speed", tag_next_).ViewComponent("cell", false);

  S_->GetEvaluator("surface-co2_concentration", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& co2 =
    *S_->Get<CompositeVector>("surface-co2_concentration", tag_next_).ViewComponent("cell", false);

  // note that this is used as the column area, which is maybe not always
  // right.  Likely correct for soil carbon calculations and incorrect for
  // surface vegetation calculations (where the subsurface's face area is more
  // correct?)
  S_->GetEvaluator("surface-cell_volume", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& scv =
    *S_->Get<CompositeVector>("surface-cell_volume", tag_next_).ViewComponent("cell", false);

  // Create workspace arrays (these should be removed when data is correctly oriented).
  auto temp_c = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto pres_c = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto dz_c = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto depth_c = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));

  // Create a workspace array for the result
  Epetra_SerialDenseVector co2_decomp_c(ncells_per_col_);
  Epetra_SerialDenseVector trans_c(ncells_per_col_);
  double sw_c(0.);

  total_lai.PutScalar(0.);

  // loop over columns and apply the model
  for (AmanziMesh::Entity_ID col = 0; col != num_cols_; ++col) {
    // update the various soil arrays
    FieldToColumn_(col, *temp(0), temp_c.ptr());
    FieldToColumn_(col, *pres(0), pres_c.ptr());
    ColDepthDz_(col, depth_c.ptr(), dz_c.ptr());

    // copy over the soil carbon arrays
    auto col_iter = mesh_->columns.getCells(col);
    ncells_per_col_ = col_iter.size();

    // -- serious cache thrash... --etc
    for (std::size_t i = 0; i != col_iter.size(); ++i) {
      AmanziGeometry::Point centroid = mesh_->getCellCentroid(col_iter[i]);
      //      std::cout << "Col iter col=" << col << ", index i=" << i << ", cell=" << col_iter[i] << " at " << centroid << std::endl;
      for (int p = 0; p != soil_carbon_pools_[col][i]->nPools; ++p) {
        soil_carbon_pools_[col][i]->SOM[p] = sc_pools[p][col_iter[i]];
      }
    }

    // Create the Met data struct
    MetData met;
    met.qSWin = qSWin[0][col];
    met.tair = air_temp[0][col];
    met.windv = wind_speed[0][col];
    met.wind_ref_ht = wind_speed_ref_ht_;
    met.vp_air = vp_air[0][col];
    met.CO2a = co2[0][col];
    met.lat = lat_;
    sw_c = met.qSWin;

    // call the model
    BGCAdvance(S_->get_time(tag_current_),
               dt,
               scv[0][col],
               cryoturbation_coef_,
               met,
               *temp_c,
               *pres_c,
               *depth_c,
               *dz_c,
               pfts_[col],
               soil_carbon_pools_[col],
               co2_decomp_c,
               trans_c,
               sw_c);

    // copy back
    // -- serious cache thrash... --etc
    for (std::size_t i = 0; i != col_iter.size(); ++i) {
      for (int p = 0; p != soil_carbon_pools_[col][i]->nPools; ++p) {
        sc_pools[p][col_iter[i]] = soil_carbon_pools_[col][i]->SOM[p];
      }

      // and integrate the decomp
      co2_decomp[0][col_iter[i]] += co2_decomp_c[i];


      // and pull in the transpiration, converting to mol/m^3/s, as a sink
      trans[0][col_iter[i]] = trans_c[i] / .01801528;
      sw[0][col] = sw_c;
    }

    for (int lcv_pft = 0; lcv_pft != pfts_[col].size(); ++lcv_pft) {
      biomass[lcv_pft][col] = pfts_[col][lcv_pft]->totalBiomass;
      leafbiomass[lcv_pft][col] = pfts_[col][lcv_pft]->Bleaf;
      csink[lcv_pft][col] = pfts_[col][lcv_pft]->CSinkLimit;
      lai[lcv_pft][col] = pfts_[col][lcv_pft]->lai;

      total_transpiration[lcv_pft][col] = pfts_[col][lcv_pft]->ET / 0.01801528;
      total_lai[0][col] += pfts_[col][lcv_pft]->lai;
    }

  } // end loop over columns

  // mark primaries as changed
  changedEvaluatorPrimary(trans_key_, tag_next_, *S_);
  changedEvaluatorPrimary(shaded_sw_key_, tag_next_, *S_);
  changedEvaluatorPrimary(total_lai_key_, tag_next_, *S_);
  return false;
}


// helper function for pushing field to column
void
BGCSimple::FieldToColumn_(AmanziMesh::Entity_ID col,
                          const Epetra_Vector& vec,
                          Teuchos::Ptr<Epetra_SerialDenseVector> col_vec,
                          bool copy)
{
  if (col_vec == Teuchos::null) {
    col_vec = Teuchos::ptr(new Epetra_SerialDenseVector(ncells_per_col_));
  }

  auto col_iter = mesh_->columns.getCells(col);
  for (std::size_t i = 0; i != col_iter.size(); ++i) { (*col_vec)[i] = vec[col_iter[i]]; }
}

// helper function for pushing field to column
void
BGCSimple::FieldToColumn_(AmanziMesh::Entity_ID col,
                          const Epetra_Vector& vec,
                          double* col_vec,
                          int ncol)
{
  auto col_iter = mesh_->columns.getCells(col);
  for (std::size_t i = 0; i != col_iter.size(); ++i) { col_vec[i] = vec[col_iter[i]]; }
}


// helper function for collecting column dz and depth
void
BGCSimple::ColDepthDz_(AmanziMesh::Entity_ID col,
                       Teuchos::Ptr<Epetra_SerialDenseVector> depth,
                       Teuchos::Ptr<Epetra_SerialDenseVector> dz)
{
  AmanziMesh::Entity_ID f_above = mesh_surf_->getEntityParent(AmanziMesh::Entity_kind::CELL, col);
  auto col_iter = mesh_->columns.getCells(col);
  ncells_per_col_ = col_iter.size();

  AmanziGeometry::Point surf_centroid = mesh_->getFaceCentroid(f_above);
  AmanziGeometry::Point neg_z(3);
  neg_z.set(0., 0., -1);

  for (std::size_t i = 0; i != col_iter.size(); ++i) {
    // depth centroid
    (*depth)[i] = surf_centroid[2] - mesh_->getCellCentroid(col_iter[i])[2];

    // dz
    // -- find face_below
    const auto& [faces, dirs] = mesh_->getCellFacesAndDirections(col_iter[i]);

    // -- mimics implementation of buildColumns() in Mesh
    double mindp = 999.0;
    AmanziMesh::Entity_ID f_below = -1;
    for (std::size_t j = 0; j != faces.size(); ++j) {
      AmanziGeometry::Point normal = mesh_->getFaceNormal(faces[j]);
      if (dirs[j] == -1) normal *= -1;
      normal /= AmanziGeometry::norm(normal);

      double dp = -normal * neg_z;
      if (dp < mindp) {
        mindp = dp;
        f_below = faces[j];
      }
    }

    // -- fill the val
    (*dz)[i] = mesh_->getFaceCentroid(f_above)[2] - mesh_->getFaceCentroid(f_below)[2];
    AMANZI_ASSERT((*dz)[i] > 0.);
    f_above = f_below;
  }
}


} // namespace BGC
} // namespace Amanzi
