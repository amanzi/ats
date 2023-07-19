/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Markus Berndt
           Daniil Svyatskiy
*/

//! Subsidence through bulk ice loss and cell volumetric change.
#include "CompositeVectorFunctionFactory.hh"
#include "pk_helpers.hh"
#include "volumetric_deformation.hh"

#define DEBUG 0

namespace Amanzi {
namespace Deform {

using namespace Amanzi::AmanziMesh;


VolumetricDeformation::VolumetricDeformation(Teuchos::ParameterList& pk_tree,
                                             const Teuchos::RCP<Teuchos::ParameterList>& glist,
                                             const Teuchos::RCP<State>& S,
                                             const Teuchos::RCP<TreeVector>& solution)
  : PK(pk_tree, glist, S, solution),
    PK_Physical_Default(pk_tree, glist, S, solution),
    surf_mesh_(Teuchos::null),
    deformed_this_step_(false)
{
  dt_max_ = plist_->get<double>("max time step [s]", std::numeric_limits<double>::max());

  // The deformation mode describes how to calculate new cell volume from a
  // provided function and the old cell volume.
  std::string mode_name = plist_->get<std::string>("deformation mode", "prescribed");
  if (mode_name == "prescribed") {
    deform_mode_ = DEFORM_MODE_DVDT;

  } else if (mode_name == "structural") {
    deform_mode_ = DEFORM_MODE_STRUCTURAL;
    deform_region_ = plist_->get<std::string>("deformation region");
    time_scale_ = plist_->get<double>("deformation relaxation time [s]", 60.);
    structural_vol_frac_ =
      plist_->get<double>("volume fraction of solids required to be structural [-]", 0.45);
    overpressured_limit_ = plist_->get<double>("overpressured relative compressibility limit", 0.2);

  } else if (mode_name == "saturation") {
    deform_mode_ = DEFORM_MODE_SATURATION;
    deform_region_ = plist_->get<std::string>("deformation region");
    min_S_liq_ = plist_->get<double>("minimum liquid saturation", 0.3);
    min_porosity_ = plist_->get<double>("minimum porosity", 0.5);
    deform_scaling_ = plist_->get<double>("deformation scaling", 1.);
    overpressured_limit_ = plist_->get<double>("overpressured relative compressibility limit", 0.2);

  } else {
    Errors::Message mesg(
      "Unknown deformation mode specified.  Valid: [prescribed, structural, saturation].");
    Exceptions::amanzi_throw(mesg);
  }

  // The deformation strategy describes how to calculate nodal deformation
  // from cell volume change.
  std::string strategy_name =
    plist_->get<std::string>("deformation strategy", "global optimization");
  if (strategy_name == "global optimization") {
    strategy_ = DEFORM_STRATEGY_GLOBAL_OPTIMIZATION;
    Errors::Message mesg("Deformation strategy \"global optimization\" is no longer supported.");
    Exceptions::amanzi_throw(mesg);
  } else if (strategy_name == "mstk implementation") {
    strategy_ = DEFORM_STRATEGY_MSTK;
  } else if (strategy_name == "average") {
    strategy_ = DEFORM_STRATEGY_AVERAGE;
  } else {
    Errors::Message mesg("Unknown deformation strategy specified. Valid: [global optimization, "
                         "mstk implementation, average]");
    Exceptions::amanzi_throw(mesg);
  }

  // collect keys
  domain_surf_ = Keys::readDomainHint(*plist_, domain_, "domain", "surface");
  domain_surf_3d_ = domain_surf_ + "_3d";

  sat_liq_key_ = Keys::readKey(*plist_, domain_, "saturation liquid", "saturation_liquid");
  sat_gas_key_ = Keys::readKey(*plist_, domain_, "saturation gas", "saturation_gas");
  sat_ice_key_ = Keys::readKey(*plist_, domain_, "saturation ice", "saturation_ice");
  cv_key_ = Keys::readKey(*plist_, domain_, "cell volume", "cell_volume");
  del_cv_key_ = Keys::readKey(*plist_, domain_, "cell volume change", "delta_cell_volume");
  poro_key_ = Keys::readKey(*plist_, domain_, "porosity", "porosity");
  vertex_loc_key_ = Keys::getKey(domain_, "vertex_coordinates");
  vertex_loc_surf3d_key_ = Keys::getKey(domain_surf_3d_, "vertex_coordinates");
  nodal_dz_key_ = Keys::readKey(*plist_, domain_, "vertex coordinate dz", "nodal_dz");
  face_above_dz_key_ = Keys::readKey(*plist_, domain_, "face above dz", "face_above_dz");
}

// -- Setup data
void
VolumetricDeformation::Setup()
{
  PK_Physical_Default::Setup();

  // Save both non-const deformable versions of the meshes and create storage
  // for the vertex coordinates.  These are saved to be able to create the
  // deformed mesh after restart, and additionally must be checkpointed.
  int dim = mesh_->getSpaceDimension();
  mesh_nc_ = S_->GetDeformableMesh(domain_);

  S_->Require<CompositeVector, CompositeVectorSpace>(vertex_loc_key_, tag_next_, vertex_loc_key_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->SetComponent("node", AmanziMesh::Entity_kind::NODE, dim);
  if (tag_next_ != Amanzi::Tags::NEXT) {
    S_->Require<CompositeVector, CompositeVectorSpace>(
        vertex_loc_key_, Amanzi::Tags::NEXT, vertex_loc_key_)
      .SetMesh(mesh_)
      ->SetGhosted()
      ->SetComponent("node", AmanziMesh::Entity_kind::NODE, dim);
  }

  // note this mesh is never deformed in this PK as all movement is vertical
  if (S_->HasMesh(domain_surf_)) surf_mesh_ = S_->GetMesh(domain_surf_);

  // the 3D variant of the surface mesh does change!
  if (S_->HasMesh(domain_surf_3d_)) {
    surf3d_mesh_ = S_->GetMesh(domain_surf_3d_);
    surf3d_mesh_nc_ = S_->GetDeformableMesh(domain_surf_3d_);
    S_->Require<CompositeVector, CompositeVectorSpace>(
        vertex_loc_surf3d_key_, tag_next_, vertex_loc_surf3d_key_)
      .SetMesh(surf3d_mesh_)
      ->SetGhosted()
      ->SetComponent("node", AmanziMesh::Entity_kind::NODE, dim);
    if (tag_next_ != Amanzi::Tags::NEXT) {
      S_->Require<CompositeVector, CompositeVectorSpace>(
          vertex_loc_surf3d_key_, Amanzi::Tags::NEXT, vertex_loc_surf3d_key_)
        .SetMesh(mesh_)
        ->SetGhosted()
        ->SetComponent("node", AmanziMesh::Entity_kind::NODE, dim);
    }
  }

  // create storage for primary variable, base_porosity
  S_->Require<CompositeVector, CompositeVectorSpace>(key_, tag_current_, name_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  // Create storage and a function for cell volume change
  //
  // Note, this is purely work space and need not be kept, but we put it in
  // state for debugging purposes.
  auto& cv_fac = S_->Require<CompositeVector, CompositeVectorSpace>(del_cv_key_, tag_next_, name_);
  cv_fac.SetMesh(mesh_)->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  switch (deform_mode_) {
  case (DEFORM_MODE_DVDT): {
    // Create the deformation function
    Teuchos::ParameterList func_plist = plist_->sublist("deformation function");
    std::vector<std::string> compnames;
    deform_func_ = Functions::CreateCompositeVectorFunction(func_plist, cv_fac, compnames);
    // note, should check that cells exist in the function?
    break;
  }

  case (DEFORM_MODE_SATURATION, DEFORM_MODE_STRUCTURAL): {
    requireAtNext(sat_liq_key_, tag_next_, *S_);
    requireAtCurrent(sat_liq_key_, tag_current_, *S_, name_)
      .SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    requireAtNext(sat_ice_key_, tag_next_, *S_);
    requireAtCurrent(sat_ice_key_, tag_current_, *S_, name_)
      .SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    requireAtNext(sat_gas_key_, tag_next_, *S_);
    requireAtCurrent(sat_gas_key_, tag_current_, *S_, name_)
      .SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    requireAtNext(poro_key_, tag_next_, *S_);
    requireAtCurrent(poro_key_, tag_current_, *S_, name_)
      .SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    break;
  }
  default: {
    AMANZI_ASSERT(0);
  }
  }

  // require for cell volume, and make sure the cell volume is deformable!
  Teuchos::ParameterList& cv_eval_list = S_->GetEvaluatorList(cv_key_);
  if (!cv_eval_list.isParameter("evaluator type")) {
    // empty list -- set it to be deforming cell volume, that uses our primary
    // variable as its indicator to change the values of cell volume.
    cv_eval_list.set<std::string>("evaluator type", "deforming cell volume");
    cv_eval_list.set<std::string>("deformation key", key_);
  }
  requireAtNext(cv_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  // require a copy at the old tag
  requireAtCurrent(cv_key_, tag_current_, *S_, name_);

  // Strategy-specific setup
  switch (strategy_) {
  case (DEFORM_STRATEGY_GLOBAL_OPTIMIZATION): {
    // // create the operator
    // Teuchos::ParameterList op_plist = plist_->sublist("global solve operator");
    // // def_matrix_ = Teuchos::rcp(new Operators::MatrixVolumetricDeformation(op_plist, mesh_));

    // // NOTE: this doesn't work because ATS operators are not supported, and
    // // this isn't based on Amanzi::Operators::Operator.  But this code branch
    // // is fairly dead anyway.... --etc
    // def_matrix_->set_inverse_parameters();

    // // create storage for the nodal deformation
    // S_->RequireField(Keys::getKey(domain_,"nodal_dz"), name_)->SetMesh(mesh_)->SetGhosted()
    //     ->SetComponent("node", AmanziMesh::Entity_kind::NODE, 1);
    break;
  }

  case (DEFORM_STRATEGY_MSTK): {
    requireAtCurrent(sat_ice_key_, tag_current_, *S_, name_)
      .SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

    requireAtCurrent(poro_key_, tag_current_, *S_, name_)
      .SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    break;
  }

  case (DEFORM_STRATEGY_AVERAGE): {
    // Both of these are work space and do not need to be put into State, but
    // we do so anyway to allow them to be visualized for debugging (hence at
    // tag_next_).

    // create storage for the nodal change in position
    //  dof 0: the volume-averaged displacement
    //  dof 1: the simply averaged displacement
    //  dof 2: count of number of faces contributing to the node change
    requireAtNext(nodal_dz_key_, tag_next_, *S_, name_)
      .SetMesh(mesh_)
      ->SetGhosted()
      ->SetComponent("node", AmanziMesh::Entity_kind::NODE, 3);

    // create cell-based storage for deformation of the face above the cell
    S_->Require<CompositeVector, CompositeVectorSpace>(face_above_dz_key_, tag_next_, name_)
      .SetMesh(mesh_)
      ->SetGhosted()
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    break;
  }
  default: {
  }
  }
}


// -- Initialize owned (dependent) variables.
void
VolumetricDeformation::Initialize()
{
  // sets base porosity
  PK_Physical_Default::Initialize();

  // initialize the deformation
  S_->GetW<CompositeVector>(del_cv_key_, tag_next_, name_).PutScalar(0.);
  S_->GetRecordW(del_cv_key_, tag_next_, name_).set_initialized();

  switch (strategy_) {
  case (DEFORM_STRATEGY_GLOBAL_OPTIMIZATION): {
    // // initialize the initial displacement to be zero
    // S_->GetW<CompositeVector>(nodal_dz_key_, tag_next_, name_).PutScalar(0.);
    // S_->GetRecordW(nodal_dz_key_, tag_next_, name_).set_initialized();
    // break;
  }
  case (DEFORM_STRATEGY_AVERAGE): {
    // initialize the initial displacement to be zero
    S_->GetW<CompositeVector>(nodal_dz_key_, tag_next_, name_).PutScalar(0.);
    S_->GetRecordW(nodal_dz_key_, tag_next_, name_).set_initialized();
    S_->GetW<CompositeVector>(face_above_dz_key_, tag_next_, name_).PutScalar(0.);
    S_->GetRecordW(face_above_dz_key_, tag_next_, name_).set_initialized();
    break;
  }
  default: {
  }
  }

  // initialize the vertex coordinate to the current mesh
  copyMeshCoordinatesToVector(
    *mesh_, S_->GetW<CompositeVector>(vertex_loc_key_, tag_next_, vertex_loc_key_));
  S_->GetRecordW(vertex_loc_key_, tag_next_, vertex_loc_key_).set_initialized();
  if (tag_next_ != Amanzi::Tags::NEXT) {
    copyMeshCoordinatesToVector(
      *mesh_, S_->GetW<CompositeVector>(vertex_loc_key_, Amanzi::Tags::NEXT, vertex_loc_key_));
    S_->GetRecordW(vertex_loc_key_, Amanzi::Tags::NEXT, vertex_loc_key_).set_initialized();
  }

  if (surf3d_mesh_ != Teuchos::null) {
    copyMeshCoordinatesToVector(
      *surf3d_mesh_,
      S_->GetW<CompositeVector>(vertex_loc_surf3d_key_, tag_next_, vertex_loc_surf3d_key_));
    S_->GetRecordW(vertex_loc_surf3d_key_, tag_next_, vertex_loc_surf3d_key_).set_initialized();
    if (tag_next_ != Amanzi::Tags::NEXT) {
      copyMeshCoordinatesToVector(*surf3d_mesh_,
                                  S_->GetW<CompositeVector>(vertex_loc_surf3d_key_,
                                                            Amanzi::Tags::NEXT,
                                                            vertex_loc_surf3d_key_));
      S_->GetRecordW(vertex_loc_surf3d_key_, Amanzi::Tags::NEXT, vertex_loc_surf3d_key_)
        .set_initialized();
    }
  }
}


bool
VolumetricDeformation::AdvanceStep(double t_old, double t_new, bool reinit)
{
  deformed_this_step_ = false;
  double dt = t_new - t_old;
  Teuchos::OSTab out = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "----------------------------------------------------------------" << std::endl
               << "Advancing: t0 = " << t_old << " t1 = " << t_new << " h = " << dt << std::endl
               << "----------------------------------------------------------------" << std::endl;

  // Collect data from state
  auto dcell_vol_vec = S_->GetPtrW<CompositeVector>(del_cv_key_, tag_next_, name_);
  dcell_vol_vec->PutScalar(0.);

  // Calculate the change in cell volumes
  switch (deform_mode_) {
  case (DEFORM_MODE_DVDT): {
    deform_func_->Compute((t_old + t_new) / 2., dcell_vol_vec.ptr());
    dcell_vol_vec->Scale(dt);
    break;
  }

  case (DEFORM_MODE_SATURATION): {
    if (S_->HasEvaluator(cv_key_, tag_current_))
      S_->GetEvaluator(cv_key_, tag_current_).Update(*S_, name_);
    if (S_->HasEvaluator(sat_liq_key_, tag_current_))
      S_->GetEvaluator(sat_liq_key_, tag_current_).Update(*S_, name_);
    if (S_->HasEvaluator(sat_gas_key_, tag_current_))
      S_->GetEvaluator(sat_gas_key_, tag_current_).Update(*S_, name_);
    if (S_->HasEvaluator(sat_ice_key_, tag_current_))
      S_->GetEvaluator(sat_ice_key_, tag_current_).Update(*S_, name_);
    if (S_->HasEvaluator(poro_key_, tag_current_))
      S_->GetEvaluator(poro_key_, tag_current_).Update(*S_, name_);

    const Epetra_MultiVector& cv =
      *S_->Get<CompositeVector>(cv_key_, tag_current_).ViewComponent("cell", true);
    const Epetra_MultiVector& s_liq =
      *S_->Get<CompositeVector>(sat_liq_key_, tag_current_).ViewComponent("cell", false);
    const Epetra_MultiVector& s_ice =
      *S_->Get<CompositeVector>(sat_ice_key_, tag_current_).ViewComponent("cell", false);
    const Epetra_MultiVector& s_gas =
      *S_->Get<CompositeVector>(sat_gas_key_, tag_current_).ViewComponent("cell", false);
    const Epetra_MultiVector& poro =
      *S_->Get<CompositeVector>(poro_key_, tag_current_).ViewComponent("cell", false);
    const Epetra_MultiVector& base_poro =
      *S_->Get<CompositeVector>(key_, tag_current_).ViewComponent("cell", false);

    Epetra_MultiVector& dcell_vol_c = *dcell_vol_vec->ViewComponent("cell", false);
    int dim = mesh_->getSpaceDimension();

    auto cells = mesh_->getSetEntities(
      deform_region_, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

    for (auto c : cells) {
      double frac = 0.;

      if (s_liq[0][c] > min_S_liq_) { // perform deformation if s_liq > min_S_liq_
        if ((poro[0][c] - base_poro[0][c]) / base_poro[0][c] <
            overpressured_limit_) { // perform deformation
          // if pressure have been relaxed enough
          frac = std::min((base_poro[0][c] - min_porosity_) / (1 - min_porosity_),
                          deform_scaling_ * ((1 - s_ice[0][c]) - min_S_liq_) * base_poro[0][c]);
        }
      }
      dcell_vol_c[0][c] = -frac * cv[0][c];

#if DEBUG
      double soil_mass_vol = cv[0][c] * (1 - base_poro[0][c]);
      std::cout << "Cell: " << c << " " << cv[0][c] << " " << dcell_vol_c[0][c] << " frac " << frac
                << " poro " << poro[0][c] << " ice " << s_ice[0][c] << " liq " << s_liq[0][c]
                << " gas " << s_gas[0][c] << " soil vol " << soil_mass_vol << std::endl;
#endif
    }
    break;
  }

  case (DEFORM_MODE_STRUCTURAL): {
    if (S_->HasEvaluator(cv_key_, tag_current_))
      S_->GetEvaluator(cv_key_, tag_current_).Update(*S_, name_);
    if (S_->HasEvaluator(sat_liq_key_, tag_current_))
      S_->GetEvaluator(sat_liq_key_, tag_current_).Update(*S_, name_);
    if (S_->HasEvaluator(sat_gas_key_, tag_current_))
      S_->GetEvaluator(sat_gas_key_, tag_current_).Update(*S_, name_);
    if (S_->HasEvaluator(sat_ice_key_, tag_current_))
      S_->GetEvaluator(sat_ice_key_, tag_current_).Update(*S_, name_);
    if (S_->HasEvaluator(poro_key_, tag_current_))
      S_->GetEvaluator(poro_key_, tag_current_).Update(*S_, name_);

    const Epetra_MultiVector& cv =
      *S_->Get<CompositeVector>(cv_key_, tag_current_).ViewComponent("cell", true);
    const Epetra_MultiVector& s_liq =
      *S_->Get<CompositeVector>(sat_liq_key_, tag_current_).ViewComponent("cell", false);
    const Epetra_MultiVector& s_ice =
      *S_->Get<CompositeVector>(sat_ice_key_, tag_current_).ViewComponent("cell", false);
    const Epetra_MultiVector& s_gas =
      *S_->Get<CompositeVector>(sat_gas_key_, tag_current_).ViewComponent("cell", false);
    const Epetra_MultiVector& poro =
      *S_->Get<CompositeVector>(poro_key_, tag_current_).ViewComponent("cell", false);
    const Epetra_MultiVector& base_poro =
      *S_->Get<CompositeVector>(key_, tag_current_).ViewComponent("cell", false);

    Epetra_MultiVector& dcell_vol_c = *dcell_vol_vec->ViewComponent("cell", false);
    dcell_vol_c.PutScalar(0.);
    int dim = mesh_->getSpaceDimension();

    auto cells = mesh_->getSetEntities(
      deform_region_, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

    double time_factor = dt > time_scale_ ? 1 : dt / time_scale_;
    for (auto c : cells) {
      double fs = 1 - poro[0][c];
      double fi = poro[0][c] * s_ice[0][c];

      double frac = 0.;
      if (fs + fi < structural_vol_frac_ && // sub-structural... start subsiding
          (poro[0][c] - base_poro[0][c]) / base_poro[0][c] < overpressured_limit_) {
        // perform deformation if we are not too overpressured
        frac = (structural_vol_frac_ - (fs + fi)) * time_factor;
      }

      dcell_vol_c[0][c] = -frac * cv[0][c];

      AMANZI_ASSERT(dcell_vol_c[0][c] <= 0);
#if DEBUG
      std::cout << "Cell " << c << ": V, dV: " << cv[0][c] << " " << dcell_vol_c[0][c] << std::endl
                << "  poro_0 " << base_poro[0][c] << " | poro " << poro[0][c] << " | frac " << frac
                << " | time factor " << time_factor << std::endl
                << "  ice " << s_ice[0][c] << " | liq " << s_liq[0][c] << " | gas " << s_gas[0][c]
                << std::endl
                << "  conserved soil vol:" << cv[0][c] * (1 - base_poro[0][c]) << std::endl;
#endif
    }
    break;
  }
  default:
    AMANZI_ASSERT(0);
  }

  // set up the fixed nodes
  Teuchos::RCP<AmanziMesh::Entity_ID_List> fixed_node_list;
  if (strategy_ == DEFORM_STRATEGY_GLOBAL_OPTIMIZATION || strategy_ == DEFORM_STRATEGY_MSTK) {
    // set up the fixed list
    fixed_node_list = Teuchos::rcp(new AmanziMesh::Entity_ID_List());
    auto nodes = mesh_->getSetEntities(
      "bottom face", AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_kind::OWNED);
    for (auto n : nodes) fixed_node_list->push_back(n);
  }

  // only deform if needed
  double dcell_vol_norm(0.);
  dcell_vol_vec->NormInf(&dcell_vol_norm);

  if (dcell_vol_norm > 0.) {
    // Deform the subsurface mesh
    switch (strategy_) {
    case (DEFORM_STRATEGY_MSTK): {
      // collect needed data, ghosted
      const Epetra_MultiVector& poro =
        *S_->Get<CompositeVector>(poro_key_, tag_current_).ViewComponent("cell");
      const Epetra_MultiVector& s_ice =
        *S_->Get<CompositeVector>(sat_ice_key_, tag_current_).ViewComponent("cell", false);
      const Epetra_MultiVector& cv =
        *S_->Get<CompositeVector>(cv_key_, tag_current_).ViewComponent("cell");

      // -- dcell vol
      const Epetra_MultiVector& dcell_vol_c = *dcell_vol_vec->ViewComponent("cell", true);

      // data needed in vectors
      int ncells = mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);
      int nnodes = mesh_->getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_kind::ALL);
      std::vector<double> target_cell_vols(ncells);
      std::vector<double> min_cell_vols(ncells);
      int dim = mesh_->getSpaceDimension();
      double min_height = std::numeric_limits<double>::max();

      for (int c = 0; c != ncells; ++c) {
        target_cell_vols[c] = cv[0][c] + dcell_vol_c[0][c];
        min_cell_vols[c] = ((1 - poro[0][c]) + poro[0][c] * s_ice[0][c]) * cv[0][c];
        // min vol is rock vol + ice + a bit
        if (std::abs(cv[0][c] - target_cell_vols[c]) / cv[0][c] > 1e-4) {
#if DEBUG
          std::cout << "Cell " << c << ": "
                    << "V, V_target, V_min " << cv[0][c] << " " << target_cell_vols[c] << " "
                    << min_cell_vols[c] << std::endl;
#endif
          auto centroid = mesh_->getCellCentroid(c);
          min_height = std::min(min_height, centroid[dim - 1]);
          AMANZI_ASSERT(min_cell_vols[c] <= target_cell_vols[c]);

        } else {
          target_cell_vols[c] = -1.; // disregard these cells
        }
      }

      // make a list of nodes below this height to keep fixed position to help MSTK
      AmanziMesh::Entity_ID_List below_node_list;
      for (unsigned int n = 0; n != nnodes; ++n) {
        AmanziGeometry::Point nc(3);
        nc = mesh_->getNodeCoordinate(n);
        if (nc[dim - 1] < min_height) below_node_list.emplace_back(n);
      }
      
      Errors::Message mesg("Volumetric Deformation not implemented in Amanzi");
      Exceptions::amanzi_throw(mesg);
      //mesh_nc_->deform(target_cell_vols, min_cell_vols, below_node_list, true);
      deformed_this_step_ = true;
      break;
    }

    case (DEFORM_STRATEGY_AVERAGE): {
      const Epetra_MultiVector& dcell_vol_c = *dcell_vol_vec->ViewComponent("cell", true);
      const Epetra_MultiVector& cv =
        *S_->Get<CompositeVector>(cv_key_, tag_current_).ViewComponent("cell");

      CompositeVector& nodal_dz_vec = S_->GetW<CompositeVector>(nodal_dz_key_, tag_next_, name_);
      { // context for vector prior to communication
        Epetra_MultiVector& nodal_dz = *nodal_dz_vec.ViewComponent("node", "true");
        nodal_dz.PutScalar(0.);

        int ncols = mesh_->columns.num_columns_owned;
        int z_index = mesh_->getSpaceDimension() - 1;
        for (int col = 0; col != ncols; ++col) {
          auto col_cells = mesh_->columns.getCells(col);
          auto col_faces = mesh_->columns.getFaces(col);
          AMANZI_ASSERT(col_faces.size() == col_cells.size() + 1);

          // iterate up the column accumulating face displacements
          double face_displacement = 0.;
          for (int ci = col_cells.size() - 1; ci >= 0; --ci) {
            int f_below = col_faces[ci + 1];
            int f_above = col_faces[ci];

            double dz =
              mesh_->getFaceCentroid(f_above)[z_index] - mesh_->getFaceCentroid(f_below)[z_index];
            face_displacement += -dz * dcell_vol_c[0][col_cells[ci]] / cv[0][col_cells[ci]];

            AMANZI_ASSERT(face_displacement >= 0.);
#if DEBUG
            if (face_displacement > 0.) {
              std::cout << "  Shifting cell " << col_cells[ci] << ", with personal displacement of "
                        << -dz * dcell_vol_c[0][col_cells[ci]] / cv[0][col_cells[ci]]
                        << " and frac " << -dcell_vol_c[0][col_cells[ci]] / cv[0][col_cells[ci]]
                        << std::endl;
            }
#endif

            // shove the face changes into the nodal averages
            auto nodes = mesh_->getFaceNodes(f_above);
            for (auto n : nodes) {
              nodal_dz[0][n] += face_displacement;
              nodal_dz[1][n] += dz;
              nodal_dz[2][n]++;
            }
          }
        }
      }

      // take the averages
      nodal_dz_vec.GatherGhostedToMaster();
      nodal_dz_vec.ScatterMasterToGhosted();

      Epetra_MultiVector& nodal_dz = *nodal_dz_vec.ViewComponent("node", "true");
      for (int n = 0; n != nodal_dz.MyLength(); ++n) {
        if (nodal_dz[2][n] > 0) {
          nodal_dz[0][n] /= nodal_dz[2][n];
          nodal_dz[1][n] /= nodal_dz[2][n];
        }
      }

      // deform the mesh
      AmanziMesh::Entity_ID_View node_ids("node_ids", nodal_dz.MyLength());
      AmanziMesh::Point_View new_positions("new_positions", nodal_dz.MyLength());
      for (int n = 0; n != nodal_dz.MyLength(); ++n) {
        node_ids[n] = n;
        new_positions[n] = mesh_->getNodeCoordinate(n);
        AMANZI_ASSERT(nodal_dz[0][n] >= 0.);
        new_positions[n][2] -= nodal_dz[0][n];
      }

      for (auto& p : new_positions) { AMANZI_ASSERT(AmanziGeometry::norm(p) >= 0.); }
      AmanziMesh::MeshAlgorithms::deform(*mesh_nc_, node_ids, new_positions);
      deformed_this_step_ = true;
      // INSERT EXTRA CODE TO UNDEFORM THE MESH FOR MIN_VOLS!
      break;
    }
    default:
      AMANZI_ASSERT(0);
    }

    // now we have to adapt the surface mesh to the new volume mesh
    // extract the correct new coordinates for the surface from the domain
    // mesh and update the surface mesh accordingly
    if (surf3d_mesh_nc_ != Teuchos::null) {
      // done on ALL to avoid lack of communication issues in deform
      int nsurfnodes =
        surf3d_mesh_nc_->getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_kind::ALL);

      AmanziMesh::Entity_ID_View surface_nodeids("surface_nodeids", nsurfnodes);
      AmanziMesh::Point_View surface_newpos("surface_newpos", nsurfnodes);

      for (int i = 0; i != nsurfnodes; ++i) {
        // get the coords of the node
        AmanziMesh::Entity_ID pnode = surf3d_mesh_nc_->getEntityParent(AmanziMesh::Entity_kind::NODE, i);
        int dim = mesh_->getSpaceDimension();
        AmanziGeometry::Point coord_domain(dim);
        coord_domain = mesh_->getNodeCoordinate(pnode);

        surface_nodeids[i] = i;
        surface_newpos[i] = coord_domain;
      }
      AmanziMesh::MeshAlgorithms::deform(*surf3d_mesh_nc_, surface_nodeids, surface_newpos);
    }

    // Note, this order is intentionally odd.  The deforming cell volume
    // evaluator uses base porosity as it's dependency.  But the reality is
    // that it does not use the value -- that comes from the mesh.  So we mark
    // the base porosity as changed, then updated the cell volume, then use
    // that to update the base porosity values itself.

    // mark base porosity has having changed
    ChangedSolutionPK(tag_next_);

    // update cell volumes
    S_->GetEvaluator(cv_key_, tag_next_).Update(*S_, name_);
    const CompositeVector& cv_vec_new = S_->Get<CompositeVector>(cv_key_, tag_next_);
    // unclear why to scatter here, maybe to pre-scatter? --ETC
    cv_vec_new.ScatterMasterToGhosted("cell");
    const Epetra_MultiVector& cv_new = *cv_vec_new.ViewComponent("cell", false);

    const Epetra_MultiVector& cv =
      *S_->Get<CompositeVector>(cv_key_, tag_current_).ViewComponent("cell", false);

    // update base porosity
    const CompositeVector& base_poro_vec_old = S_->Get<CompositeVector>(key_, tag_current_);
    const Epetra_MultiVector& base_poro_old = *base_poro_vec_old.ViewComponent("cell", false);
    Epetra_MultiVector& base_poro =
      *S_->GetW<CompositeVector>(key_, tag_next_, name_).ViewComponent("cell");

    int ncells = base_poro.MyLength();
    for (int c = 0; c != ncells; ++c) {
      base_poro[0][c] = 1. - (1. - base_poro_old[0][c]) * cv[0][c] / cv_new[0][c];
#if DEBUG
      if (std::abs(cv_new[0][c] - cv[0][c]) > 1.e-12) {
        std::cout << "Deformed Cell " << c << ": V,V_new " << cv[0][c] << " " << cv_new[0][c]
                  << std::endl
                  << "             result porosity " << base_poro_old[0][c] << " "
                  << base_poro[0][c] << std::endl;
      }
#endif
    }

  } // if any cell volumes have changed

  // debug...
  std::vector<std::string> names = {
    "base_poro old", "base_poro new", "cell_vol old", "cell_vol new"
  };
  std::vector<Teuchos::Ptr<const CompositeVector>> vecs = {
    S_->GetPtr<CompositeVector>(key_, tag_current_).ptr(),
    S_->GetPtr<CompositeVector>(key_, tag_next_).ptr(),
    S_->GetPtr<CompositeVector>(cv_key_, tag_current_).ptr(),
    S_->GetPtr<CompositeVector>(cv_key_, tag_next_).ptr()
  };
  db_->WriteVectors(names, vecs);

  return false;
}


void
VolumetricDeformation::CommitStep(double t_old, double t_new, const Tag& tag_next)
{
  // saves primary variable
  PK_Physical_Default::CommitStep(t_old, t_new, tag_next);

  AMANZI_ASSERT(tag_next == tag_next_ || tag_next == Tags::NEXT);
  Tag tag_current = tag_next == tag_next_ ? tag_current_ : Tags::CURRENT;

  // also save conserved quantity and saturation
  if (deform_mode_ == DEFORM_MODE_SATURATION || deform_mode_ == DEFORM_MODE_STRUCTURAL) {
    assign(sat_liq_key_, tag_current, tag_next, *S_);
    assign(sat_ice_key_, tag_current, tag_next, *S_);
    assign(sat_gas_key_, tag_current, tag_next, *S_);
    assign(poro_key_, tag_current, tag_next, *S_);
  }

  if (strategy_ == DEFORM_STRATEGY_MSTK) { assign(poro_key_, tag_current, tag_next, *S_); }

  // lastly, save the new coordinates for checkpointing
  if (deformed_this_step_) {
    copyMeshCoordinatesToVector(
      *mesh_, S_->GetW<CompositeVector>(vertex_loc_key_, tag_next, vertex_loc_key_));
    if (surf3d_mesh_ != Teuchos::null) {
      copyMeshCoordinatesToVector(
        *surf3d_mesh_,
        S_->GetW<CompositeVector>(vertex_loc_surf3d_key_, tag_next, vertex_loc_surf3d_key_));
    }
  }
  assign(cv_key_, tag_current, tag_next, *S_);
}

void
VolumetricDeformation::FailStep(double t_old, double t_new, const Tag& tag)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) *vo_->os() << "Failing step." << std::endl;

  if (deformed_this_step_) {
    copyVectorToMeshCoordinates(S_->Get<CompositeVector>(vertex_loc_key_, tag), *mesh_nc_);
    if (surf3d_mesh_ != Teuchos::null) {
      copyVectorToMeshCoordinates(S_->Get<CompositeVector>(vertex_loc_surf3d_key_, tag),
                                  *surf3d_mesh_nc_);
    }
  }
}

} // namespace Deform
} // namespace Amanzi
