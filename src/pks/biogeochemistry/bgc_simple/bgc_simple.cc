/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon, Chonggang Xu

   Simple implementation of CLM's Century model for carbon decomposition and a
   simplified 2-PFT (sedge, moss) vegetation model for creating carbon.

   CURRENT ASSUMPTIONS:
     1. parallel decomp not in the vertical
     2. fields are not ordered along the column, and so must be copied
     3. all columns have the same number of cells
   ------------------------------------------------------------------------- */

#include "MeshPartition.hh"

#include "bgc_simple_funcs.hh"

#include "bgc_simple.hh"

namespace Amanzi {
namespace BGC {


BGCSimple::BGCSimple(Teuchos::ParameterList& pk_tree,
                     const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                     const Teuchos::RCP<State>& S,
                     const Teuchos::RCP<TreeVector>& solution):
  PK_Physical_Default(pk_tree, global_list, S, solution),
  PK(pk_tree, global_list, S, solution),
  ncells_per_col_(-1) {

  // set up additional primary variables -- this is very hacky...
  // -- surface energy source
  domain_surf_ = plist_->get<std::string>("surface domain name", "surface");

  // set up additional primary variables -- this is very hacky...
  Teuchos::ParameterList& FElist = S->FEList();
  // -- transpiration
  trans_key_ = Keys::readKey(*plist_, domain_, "transpiration", "transpiration");
  Teuchos::ParameterList& trans_sublist =
      FElist.sublist(trans_key_);
  trans_sublist.set("field evaluator type", "primary variable");

  // -- shortwave incoming shading
  shaded_sw_key_ = Keys::readKey(*plist_, domain_surf_, "shaded shortwave radiation", "shaded_shortwave_radiation");
  Teuchos::ParameterList& sw_sublist =
      FElist.sublist(shaded_sw_key_);
  sw_sublist.set("field evaluator type", "primary variable");

  // -- lai
  total_lai_key_ = Keys::readKey(*plist_, domain_surf_, "total leaf area index", "total_leaf_area_index");
  Teuchos::ParameterList& lai_sublist =
      FElist.sublist(total_lai_key_);
  lai_sublist.set("field evaluator type", "primary variable");
}

// is a PK
// -- Setup data
void BGCSimple::Setup(const Teuchos::Ptr<State>& S) {
  PK_Physical_Default::Setup(S);

  // initial timestep
  dt_ = plist_->get<double>("initial time step", 1.);

  // my mesh is the subsurface mesh, but we need the surface mesh, index by column, as well
  mesh_surf_ = S->GetMesh("surface");

  // Create the additional, non-managed data structures
  int nPools = plist_->get<int>("number of carbon pools", 7);

  // -- SoilCarbonParameters
  Teuchos::ParameterList& sc_params = plist_->sublist("soil carbon parameters");
  std::string mesh_part_name = sc_params.get<std::string>("mesh partition");
  const Functions::MeshPartition& mp = *S->GetMeshPartition(mesh_part_name);
  const std::vector<std::string>& regions = mp.regions();

  for (std::vector<std::string>::const_iterator region=regions.begin();
       region!=regions.end(); ++region) {
    sc_params_.push_back(Teuchos::rcp(
        new SoilCarbonParameters(nPools, sc_params.sublist(*region))));
  }

  // -- PFTs -- old and new!
  Teuchos::ParameterList& pft_params = plist_->sublist("pft parameters");
  std::vector<std::string> pft_names;
  for (Teuchos::ParameterList::ConstIterator lcv=pft_params.begin();
       lcv!=pft_params.end(); ++lcv) {
    std::string pft_name = lcv->first;
    pft_names.push_back(pft_name);
  }

  int ncols = mesh_surf_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  pfts_old_.resize(ncols);
  pfts_.resize(ncols);
  for (unsigned int col=0; col!=ncols; ++col) {
    int f = mesh_surf_->entity_get_parent(AmanziMesh::CELL, col);
    ColIterator col_iter(*mesh_, f);
    std::size_t ncol_cells = col_iter.size();

    // unclear which this should be:
    // -- col area is the true face area
    double col_area = mesh_->face_area(f);
    // -- col area is the projected face area
    // double col_area = mesh_surf_->cell_volume(col);

    if (ncells_per_col_ < 0) {
      ncells_per_col_ = ncol_cells;
    } else {
      AMANZI_ASSERT(ncol_cells == ncells_per_col_);
    }

    pfts_old_[col].resize(pft_names.size());
    pfts_[col].resize(pft_names.size());

    for (int i=0; i!=pft_names.size(); ++i) {
      std::string pft_name = pft_names[i];
      Teuchos::ParameterList& pft_plist = pft_params.sublist(pft_name);
      pfts_old_[col][i] = Teuchos::rcp(new PFT(pft_name, ncol_cells));
      pfts_old_[col][i]->Init(pft_plist,col_area);
      pfts_[col][i] = Teuchos::rcp(new PFT(*pfts_old_[col][i]));
    }
  }

  // -- soil carbon pools
  soil_carbon_pools_.resize(ncols);
  for (unsigned int col=0; col!=ncols; ++col) {
    soil_carbon_pools_[col].resize(ncells_per_col_);
    ColIterator col_iter(*mesh_, mesh_surf_->entity_get_parent(AmanziMesh::CELL, col), ncells_per_col_);

    for (std::size_t i=0; i!=col_iter.size(); ++i) {
      // col_iter[i] = cell id, mp[cell_id] = index into partition list, sc_params_[index] = correct params
      soil_carbon_pools_[col][i] = Teuchos::rcp(new SoilCarbon(sc_params_[mp[col_iter[i]]]));
    }
  }

  // requirements: primary variable
  S->RequireField(key_, name_)->SetMesh(mesh_)
      ->SetComponent("cell", AmanziMesh::CELL, nPools);

  // requirements: other primary variables
  S->RequireField(trans_key_, name_)->SetMesh(mesh_)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator(trans_key_);
  trans_eval_ = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(
      S->GetFieldEvaluator(trans_key_));
  if (trans_eval_ == Teuchos::null) {
    Errors::Message message("BGC: error, failure to initialize primary variable for transpiration");
    Exceptions::amanzi_throw(message);
  }

  S->RequireField(shaded_sw_key_, name_)->SetMesh(mesh_surf_)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator(shaded_sw_key_);
  sw_eval_ = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(
      S->GetFieldEvaluator(shaded_sw_key_));
  if (sw_eval_ == Teuchos::null) {
    Errors::Message message("BGC: error, failure to initialize primary variable for shaded shortwave");
    Exceptions::amanzi_throw(message);
  }

  S->RequireField(total_lai_key_, name_)->SetMesh(mesh_surf_)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator(total_lai_key_);
  lai_eval_ = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(
      S->GetFieldEvaluator(total_lai_key_));
  if (lai_eval_ == Teuchos::null) {
    Errors::Message message("BGC: error, failure to initialize primary variable for LAI");
    Exceptions::amanzi_throw(message);
  }
  
  // requirement: diagnostics
  S->RequireField("co2_decomposition", name_)->SetMesh(mesh_)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->RequireField("surface-total_biomass", name_)->SetMesh(mesh_surf_)
      ->SetComponent("cell", AmanziMesh::CELL, pft_names.size());
  S->RequireField("surface-leaf_biomass", name_)->SetMesh(mesh_surf_)
      ->SetComponent("cell", AmanziMesh::CELL, pft_names.size());
  S->RequireField("surface-leaf_area_index", name_)->SetMesh(mesh_surf_)
      ->SetComponent("cell", AmanziMesh::CELL, pft_names.size());
  S->RequireField("surface-c_sink_limit", name_)->SetMesh(mesh_surf_)
      ->SetComponent("cell", AmanziMesh::CELL, pft_names.size());
  S->RequireField("surface-veg_total_transpiration", name_)->SetMesh(mesh_surf_)
      ->SetComponent("cell", AmanziMesh::CELL, pft_names.size());

  // requirement: temp of each cell
  S->RequireFieldEvaluator("temperature");
  S->RequireField("temperature")->SetMesh(mesh_)->AddComponent("cell", AmanziMesh::CELL, 1);

  // requirement: pressure
  S->RequireFieldEvaluator("pressure");
  S->RequireField("pressure")->SetMesh(mesh_)->AddComponent("cell", AmanziMesh::CELL, 1);

  // requirements: surface cell volume
  S->RequireField("surface-cell_volume")->SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator("surface-cell_volume");

  // requirements: Met data
  S->RequireFieldEvaluator("surface-incoming_shortwave_radiation");
  S->RequireField("surface-incoming_shortwave_radiation")->SetMesh(mesh_surf_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator("surface-air_temperature");
  S->RequireField("surface-air_temperature")->SetMesh(mesh_surf_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator("surface-relative_humidity");
  S->RequireField("surface-relative_humidity")->SetMesh(mesh_surf_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator("surface-wind_speed");
  S->RequireField("surface-wind_speed")->SetMesh(mesh_surf_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator("surface-co2_concentration");
  S->RequireField("surface-co2_concentration")->SetMesh(mesh_surf_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  // parameters
  lat_ = plist_->get<double>("latitude [degrees]", 60.);
  wind_speed_ref_ht_ = plist_->get<double>("wind speed reference height [m]", 2.0);
  cryoturbation_coef_ = plist_->get<double>("cryoturbation mixing coefficient [cm^2/yr]", 5.0);
  cryoturbation_coef_ /= 365.25e4; // convert to m^2/day

}

// -- Initialize owned (dependent) variables.
void BGCSimple::Initialize(const Teuchos::Ptr<State>& S) {
  PK_Physical_Default::Initialize(S);

  // diagnostic variable
  S->GetFieldData("co2_decomposition", name_)->PutScalar(0.);
  S->GetField("co2_decomposition", name_)->set_initialized();

  S->GetFieldData("surface-c_sink_limit", name_)->PutScalar(0.);
  S->GetField("surface-c_sink_limit", name_)->set_initialized();

  S->GetFieldData(trans_key_, name_)->PutScalar(0.);
  S->GetField(trans_key_, name_)->set_initialized();

  S->GetFieldData("surface-total_biomass", name_)->PutScalar(0.);
  S->GetField("surface-total_biomass", name_)->set_initialized();

  S->GetFieldData("surface-leaf_area_index", name_)->PutScalar(0.);
  S->GetField("surface-leaf_area_index", name_)->set_initialized();

  S->GetFieldData("surface-veg_total_transpiration", name_)->PutScalar(0.);
  S->GetField("surface-veg_total_transpiration", name_)->set_initialized();
  
  // potentially initial aboveground vegetation data
  Teuchos::RCP<Field> leaf_biomass_field = S->GetField("surface-leaf_biomass", name_);

  // -- set the subfield names
  Teuchos::RCP<Field_CompositeVector> leaf_biomass_field_cv =
      Teuchos::rcp_dynamic_cast<Field_CompositeVector>(leaf_biomass_field);
  AMANZI_ASSERT(leaf_biomass_field_cv != Teuchos::null);

  int npft = pfts_old_[0].size();
  std::vector<std::vector<std::string> > names;
  names.resize(1);
  names[0].resize(npft);
  for (int i=0; i!=npft; ++i) names[0][i] = pfts_old_[0][i]->pft_type;
  leaf_biomass_field_cv->set_subfield_names(names);

  if (!leaf_biomass_field->initialized()) {
    // -- Calculate the IC.
    if (plist_->isSublist("leaf biomass initial condition")) {
      Teuchos::ParameterList ic_plist = plist_->sublist("leaf biomass initial condition");
      leaf_biomass_field->Initialize(ic_plist);

      // -- copy into PFTs
      Epetra_MultiVector& bio = *S->GetFieldData("surface-leaf_biomass", name_)
          ->ViewComponent("cell", false);
      
      int ncols = mesh_surf_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
      for (int col=0; col!=ncols; ++col) {
        for (int i=0; i!=npft; ++i) {
          pfts_old_[col][i]->Bleaf = bio[i][col];
        }
      }
      leaf_biomass_field->set_initialized();
    }
    
    if (!leaf_biomass_field->initialized()) {
      S->GetFieldData("surface-leaf_biomass", name_)->PutScalar(0.);
      leaf_biomass_field->set_initialized();
    }
  }
  
  // init root carbon
  Teuchos::RCP<Epetra_SerialDenseVector> col_temp =
      Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  Teuchos::RCP<Epetra_SerialDenseVector> col_depth =
      Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  Teuchos::RCP<Epetra_SerialDenseVector> col_dz =
      Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));

  S->GetFieldEvaluator("temperature")->HasFieldChanged(S, name_);
  const Epetra_Vector& temp = *(*S->GetFieldData("temperature")
				->ViewComponent("cell",false))(0);

  int ncols = mesh_surf_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  for (int col=0; col!=ncols; ++col) {
    FieldToColumn_(col, temp, col_temp.ptr());
    ColDepthDz_(col, col_depth.ptr(), col_dz.ptr());

    for (int i=0; i!=npft; ++i) {
      pfts_old_[col][i]->InitRoots(*col_temp, *col_depth, *col_dz);
    }
  }

  // ensure all initialization in both PFTs?  Not sure this is
  // necessary -- likely done in initial call to commit-state --etc
  for (int col=0; col!=ncols; ++col) {
    for (int i=0; i!=npft; ++i) {
      *pfts_[col][i] = *pfts_old_[col][i];
    }
  }
}

  
// -- Commit any secondary (dependent) variables.
void BGCSimple::CommitStep(double told, double tnew, const Teuchos::RCP<State>& S) {
  // Copy the PFT over, which includes all additional state required, commit
  // the step as succesful.
  double dt = tnew - told;

  int ncols = mesh_surf_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int npft = pfts_old_[0].size();
  for (int col=0; col!=ncols; ++col) {
    for (int i=0; i!=npft; ++i) {
      *pfts_old_[col][i] = *pfts_[col][i];
    }
  }
}

// -- advance the model
bool BGCSimple::AdvanceStep(double t_old, double t_new, bool reinit) {

  double dt = t_new - t_old;

  Teuchos::OSTab out = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "----------------------------------------------------------------" << std::endl
               << "Advancing: t0 = " << S_inter_->time()
               << " t1 = " << S_next_->time() << " h = " << dt << std::endl
               << "----------------------------------------------------------------" << std::endl;

  // Copy the PFT from old to new, in case we failed the previous attempt at
  // this timestep.  This is hackery to get around the fact that PFTs are not
  // (but should be) in state.
  AmanziMesh::Entity_ID ncols = mesh_surf_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  for (AmanziMesh::Entity_ID col=0; col!=ncols; ++col) {
    int npft = pfts_old_[col].size();
    for (int i=0; i!=npft; ++i) {
      *pfts_[col][i] = *pfts_old_[col][i];
    }
  }

  // grab the required fields
  Epetra_MultiVector& sc_pools = *S_next_->GetFieldData(key_, name_)
      ->ViewComponent("cell",false);
  Epetra_MultiVector& co2_decomp = *S_next_->GetFieldData("co2_decomposition", name_)
      ->ViewComponent("cell",false);
  Epetra_MultiVector& trans = *S_next_->GetFieldData(trans_key_, name_)
      ->ViewComponent("cell",false);
  Epetra_MultiVector& sw = *S_next_->GetFieldData(shaded_sw_key_, name_)
      ->ViewComponent("cell",false);
  Epetra_MultiVector& biomass = *S_next_->GetFieldData("surface-total_biomass", name_)
      ->ViewComponent("cell",false);
  Epetra_MultiVector& leafbiomass = *S_next_->GetFieldData("surface-leaf_biomass", name_)
      ->ViewComponent("cell",false);
  Epetra_MultiVector& csink = *S_next_->GetFieldData("surface-c_sink_limit", name_)
      ->ViewComponent("cell",false);
  Epetra_MultiVector& total_lai = *S_next_->GetFieldData(total_lai_key_, name_)
      ->ViewComponent("cell",false);
  Epetra_MultiVector& lai = *S_next_->GetFieldData("surface-leaf_area_index", name_)
      ->ViewComponent("cell",false);
  Epetra_MultiVector& veg_total_transpiration = *S_next_->GetFieldData("surface-veg_total_transpiration", name_)
      ->ViewComponent("cell",false);

  S_next_->GetFieldEvaluator("temperature")->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& temp = *S_inter_->GetFieldData("temperature")
      ->ViewComponent("cell",false);

  S_next_->GetFieldEvaluator("pressure")->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& pres = *S_inter_->GetFieldData("pressure")
      ->ViewComponent("cell",false);

  S_next_->GetFieldEvaluator("surface-incoming_shortwave_radiation")->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& qSWin = *S_next_->GetFieldData("surface-incoming_shortwave_radiation")
      ->ViewComponent("cell",false);

  S_next_->GetFieldEvaluator("surface-air_temperature")->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& air_temp = *S_next_->GetFieldData("surface-air_temperature")
      ->ViewComponent("cell",false);

  S_next_->GetFieldEvaluator("surface-relative_humidity")->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& rel_hum = *S_next_->GetFieldData("surface-relative_humidity")
      ->ViewComponent("cell",false);

  S_next_->GetFieldEvaluator("surface-wind_speed")->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& wind_speed = *S_next_->GetFieldData("surface-wind_speed")
      ->ViewComponent("cell",false);

  S_next_->GetFieldEvaluator("surface-co2_concentration")->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& co2 = *S_next_->GetFieldData("surface-co2_concentration")
      ->ViewComponent("cell",false);

  // note that this is used as the column area, which is maybe not always
  // right.  Likely correct for soil carbon calculations and incorrect for
  // surface vegetation calculations (where the subsurface's face area is more
  // correct?)
  S_inter_->GetFieldEvaluator("surface-cell_volume")->HasFieldChanged(S_inter_.ptr(), name_);
  const Epetra_MultiVector& scv = *S_inter_->GetFieldData("surface-cell_volume")
      ->ViewComponent("cell", false);

  // Create workspace arrays (these should be removed when data is correctly oriented).
  Teuchos::RCP<Epetra_SerialDenseVector> temp_c =
      Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  Teuchos::RCP<Epetra_SerialDenseVector> pres_c =
      Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  Teuchos::RCP<Epetra_SerialDenseVector> dz_c =
      Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  Teuchos::RCP<Epetra_SerialDenseVector> depth_c =
      Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));

  // Create a workspace array for the result
  Epetra_SerialDenseVector co2_decomp_c(ncells_per_col_);
  Epetra_SerialDenseVector trans_c(ncells_per_col_);
  double sw_c(0.);

  total_lai.PutScalar(0.);

  // loop over columns and apply the model
  for (AmanziMesh::Entity_ID col=0; col!=ncols; ++col) {
    // update the various soil arrays
    FieldToColumn_(col, *temp(0), temp_c.ptr());
    FieldToColumn_(col, *pres(0), pres_c.ptr());
    ColDepthDz_(col, depth_c.ptr(), dz_c.ptr());

    // copy over the soil carbon arrays
    ColIterator col_iter(*mesh_, mesh_surf_->entity_get_parent(AmanziMesh::CELL, col), ncells_per_col_);
    // -- serious cache thrash... --etc
    for (std::size_t i=0; i!=col_iter.size(); ++i) {
      AmanziGeometry::Point centroid = mesh_->cell_centroid(col_iter[i]);
	  //      std::cout << "Col iter col=" << col << ", index i=" << i << ", cell=" << col_iter[i] << " at " << centroid << std::endl;
      for (int p=0; p!=soil_carbon_pools_[col][i]->nPools; ++p) {
        soil_carbon_pools_[col][i]->SOM[p] = sc_pools[p][col_iter[i]];
      }
    }

    // Create the Met data struct
    MetData met;
    met.qSWin = qSWin[0][col];
    met.tair = air_temp[0][col];
    met.windv = wind_speed[0][col];
    met.wind_ref_ht = wind_speed_ref_ht_;
    met.relhum = rel_hum[0][col];
    met.CO2a = co2[0][col];
    met.lat = lat_;
    sw_c = met.qSWin;

    // call the model
    BGCAdvance(S_inter_->time(), dt, scv[0][col], cryoturbation_coef_, met,
               *temp_c, *pres_c, *depth_c, *dz_c,
               pfts_[col], soil_carbon_pools_[col],
               co2_decomp_c, trans_c, sw_c);

    // copy back
    // -- serious cache thrash... --etc
    for (std::size_t i=0; i!=col_iter.size(); ++i) {
      for (int p=0; p!=soil_carbon_pools_[col][i]->nPools; ++p) {
        sc_pools[p][col_iter[i]] = soil_carbon_pools_[col][i]->SOM[p];
      }

      // and integrate the decomp
      co2_decomp[0][col_iter[i]] += co2_decomp_c[i];


      // and pull in the transpiration, converting to mol/m^3/s, as a sink
      trans[0][col_iter[i]] = -trans_c[i]/ .01801528;
      //      std::cout << std::scientific;
      //      std::cout << "Transpiration at " << col_iter[i] << "," << i << " = " << trans[0][col_iter[i]] << std::endl;
      sw[0][col] = sw_c;
    }

    for (int lcv_pft=0; lcv_pft!=pfts_[col].size(); ++lcv_pft) {
      biomass[lcv_pft][col] = pfts_[col][lcv_pft]->totalBiomass;
      leafbiomass[lcv_pft][col] = pfts_[col][lcv_pft]->Bleaf;
      csink[lcv_pft][col] = pfts_[col][lcv_pft]->CSinkLimit;
      lai[lcv_pft][col] = pfts_[col][lcv_pft]->lai;

      veg_total_transpiration[lcv_pft][col] = pfts_[col][lcv_pft]->ET / 0.01801528;

      total_lai[0][col] += pfts_[col][lcv_pft]->lai;
    }

  } // end loop over columns

  // mark primaries as changed
  trans_eval_->SetFieldAsChanged(S_next_.ptr());
  sw_eval_->SetFieldAsChanged(S_next_.ptr());
  lai_eval_->SetFieldAsChanged(S_next_.ptr());
  return false;
}


// helper function for pushing field to column
void BGCSimple::FieldToColumn_(AmanziMesh::Entity_ID col, const Epetra_Vector& vec,
        Teuchos::Ptr<Epetra_SerialDenseVector> col_vec, bool copy) {
  if (col_vec == Teuchos::null) {
    col_vec = Teuchos::ptr(new Epetra_SerialDenseVector(ncells_per_col_));
  }

  ColIterator col_iter(*mesh_, mesh_surf_->entity_get_parent(AmanziMesh::CELL, col), ncells_per_col_);
  for (std::size_t i=0; i!=col_iter.size(); ++i) {
    (*col_vec)[i] = vec[col_iter[i]];
  }
}

// helper function for pushing field to column
void BGCSimple::FieldToColumn_(AmanziMesh::Entity_ID col, const Epetra_Vector& vec,
                               double* col_vec, int ncol) {
  ColIterator col_iter(*mesh_, mesh_surf_->entity_get_parent(AmanziMesh::CELL, col), ncol);
  for (std::size_t i=0; i!=col_iter.size(); ++i) {
    col_vec[i] = vec[col_iter[i]];
  }
}

  
// helper function for collecting column dz and depth
void BGCSimple::ColDepthDz_(AmanziMesh::Entity_ID col,
                            Teuchos::Ptr<Epetra_SerialDenseVector> depth,
                            Teuchos::Ptr<Epetra_SerialDenseVector> dz) {
  AmanziMesh::Entity_ID f_above = mesh_surf_->entity_get_parent(AmanziMesh::CELL, col);
  ColIterator col_iter(*mesh_, f_above, ncells_per_col_);

  AmanziGeometry::Point surf_centroid = mesh_->face_centroid(f_above);
  AmanziGeometry::Point neg_z(3);
  neg_z.set(0.,0.,-1);

  for (std::size_t i=0; i!=col_iter.size(); ++i) {
    // depth centroid
    (*depth)[i] = surf_centroid[2] - mesh_->cell_centroid(col_iter[i])[2];

    // dz
    // -- find face_below
    AmanziMesh::Entity_ID_List faces;
    std::vector<int> dirs;
    mesh_->cell_get_faces_and_dirs(col_iter[i], &faces, &dirs);

    // -- mimics implementation of build_columns() in Mesh
    double mindp = 999.0;
    AmanziMesh::Entity_ID f_below = -1;
    for (std::size_t j=0; j!=faces.size(); ++j) {
      AmanziGeometry::Point normal = mesh_->face_normal(faces[j]);
      if (dirs[j] == -1) normal *= -1;
      normal /= AmanziGeometry::norm(normal);

      double dp = -normal * neg_z;
      if (dp < mindp) {
        mindp = dp;
        f_below = faces[j];
      }
    }

    // -- fill the val
    (*dz)[i] = mesh_->face_centroid(f_above)[2] - mesh_->face_centroid(f_below)[2];
    AMANZI_ASSERT( (*dz)[i] > 0. );
    f_above = f_below;
  }

}




} // namespace
} // namespace
