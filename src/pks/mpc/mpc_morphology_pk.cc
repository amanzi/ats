/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  This is the mpc_pk component of the Amanzi code.

*/

#include "mpc_morphology_pk.hh"
#include "Mesh.hh"

namespace Amanzi {

Morphology_PK::Morphology_PK(Teuchos::ParameterList& pk_tree_or_fe_list,
                             const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                             const Teuchos::RCP<State>& S,
                             const Teuchos::RCP<TreeVector>& soln)
  : PK(pk_tree_or_fe_list, global_list, S, soln),
    PK_MPCSubcycled_ATS(pk_tree_or_fe_list, global_list, S, soln)
{
  // Create verbosity object.
  vo_ = Teuchos::null;
  Teuchos::ParameterList vlist;
  vlist.sublist("verbose object") = plist_->sublist("verbose object");
  vo_ = Teuchos::rcp(new VerboseObject("Morphology_PK", vlist));
  domain_ = plist_->get<std::string>("domain name", "domain");
  name_ = "morphology pk";

  Teuchos::Array<std::string> pk_order = plist_->get<Teuchos::Array<std::string>>("PKs order");
}

// -----------------------------------------------------------------------------
// Calculate the min of sub PKs timestep sizes.
// -----------------------------------------------------------------------------
double
Morphology_PK::get_dt()
{
  if (dt_MPC_ < 0) {
    double dt = Amanzi::PK_MPCSubcycled_ATS::get_dt();
    set_dt(dt);
    return dt;
  } else {
    return dt_MPC_;
  }
}


void
Morphology_PK::Setup(const Teuchos::Ptr<State>& S)
{
  //passwd_ = "coupled_transport";  // owner's password
  //passwd_ = "state";  // owner's password


  dt_MPC_ = plist_->get<double>("dt MPC", 31557600);
  MSF_ = plist_->get<double>("morphological scaling factor", 1);

  Amanzi::PK_MPCSubcycled_ATS::Setup(S);

  mesh_ = S->GetDeformableMesh(domain_);
  vertex_coord_key_ = Keys::getKey(domain_, "vertex_coordinate");
  if (domain_ == "surface") {
    domain_3d_ = "surface_3d";
    domain_ss_ = "domain";
    vertex_coord_key_3d_ = Keys::getKey(domain_3d_, "vertex_coordinate");
    vertex_coord_key_ss_ = Keys::getKey(domain_ss_, "vertex_coordinate");
    mesh_3d_ = S->GetDeformableMesh(domain_ + "_3d");
    mesh_ss_ = S->GetDeformableMesh(domain_ss_);
  }

  // create storage for the vertex coordinates
  // we need to checkpoint those to be able to create
  // the deformed mesh after restart
  std::vector<AmanziMesh::Entity_kind> location(1);
  std::vector<int> num_dofs(1);
  std::vector<std::string> name(1);

  if (!S->HasField(vertex_coord_key_)) {
    int dim = mesh_->getSpaceDimension();
    location[0] = AmanziMesh::Entity_kind::NODE;
    num_dofs[0] = dim;
    name[0] = "node";

    S->Require<CompositeVector, CompositeVectorSpace>(vertex_coord_key_, Tags::NEXT, "state")
      .SetMesh(mesh_)
      ->SetGhosted()
      ->SetComponents(name, location, num_dofs);
  }

  if (S->HasMesh(domain_3d_) && (!S->HasField(vertex_coord_key_3d_))) {
    int dim = mesh_3d_->getSpaceDimension();
    location[0] = AmanziMesh::Entity_kind::NODE;
    num_dofs[0] = dim;
    name[0] = "node";

    S->Require<CompositeVector, CompositeVectorSpace>(vertex_coord_key_3d_, Tags::NEXT, "state")
      .SetMesh(mesh_3d_)
      ->SetGhosted()
      ->SetComponents(name, location, num_dofs);
  }

  if (S->HasMesh(domain_ss_) && (!S->HasField(vertex_coord_key_ss_))) {
    int dim = mesh_ss_->getSpaceDimension();
    location[0] = AmanziMesh::Entity_kind::NODE;
    num_dofs[0] = dim;
    name[0] = "node";

    S->Require<CompositeVector, CompositeVectorSpace>(vertex_coord_key_ss_, Tags::NEXT, "state")
      .SetMesh(mesh_ss_)
      ->SetGhosted()
      ->SetComponents(name, location, num_dofs);
  }

  elevation_increase_key_ = Keys::getKey(domain_, "deformation");
  if (!S->HasField(elevation_increase_key_)) {
    S->Require<CompositeVector, CompositeVectorSpace>(elevation_increase_key_, Tags::NEXT, "state")
      .SetMesh(mesh_)
      ->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    Teuchos::ParameterList deform_plist;
    deform_plist.set("evaluator name", elevation_increase_key_);
    deform_eval_ = Teuchos::rcp(new EvaluatorPrimary(deform_plist));
    S->SetEvaluator(elevation_increase_key_, deform_eval_);
  }

  int num_veg_species = plist_->get<int>("number of vegitation species", 1);


  Key biomass_key = Keys::getKey(domain_, "biomass");
  Key stem_density_key = Keys::getKey(domain_, "stem_density");
  Key stem_height_key = Keys::getKey(domain_, "stem_height");
  Key stem_diameter_key = Keys::getKey(domain_, "stem_diameter");
  Key plant_area_key = Keys::getKey(domain_, "plant_area");

  location[0] = AmanziMesh::Entity_kind::CELL;
  num_dofs[0] = num_veg_species;
  name[0] = "cell";

  if (!S->HasField(biomass_key)) {
    S->Require<CompositeVector, CompositeVectorSpace>(biomass_key, Tags::NEXT, biomass_key)
      .SetMesh(mesh_)
      ->SetGhosted()
      ->SetComponents(name, location, num_dofs);
    S->RequireEvaluator(biomass_key);
  }

  if (!S->HasField("msl")) S->RequireScalar("msl");
}

void
Morphology_PK::Initialize(const Teuchos::Ptr<State>& S)
{
  Amanzi::PK_MPCSubcycled_ATS::Initialize(S);

  // initialize the vertex coordinate of existing meshes

  if (S->HasField(vertex_coord_key_)) Initialize_MeshVertices_(S, mesh_, vertex_coord_key_);

  if (S->HasField(vertex_coord_key_3d_))
    Initialize_MeshVertices_(S, mesh_3d_, vertex_coord_key_3d_);

  if (S->HasField(vertex_coord_key_ss_))
    Initialize_MeshVertices_(S, mesh_ss_, vertex_coord_key_ss_);

  if (S->HasField(elevation_increase_key_)) {
    S->GetW<CompositeVector>(elevation_increase_key_, "state").PutScalar(0.);
    S->GetField(elevation_increase_key_, "state")->set_initialized();
  }

  flow_pk_ = Teuchos::rcp_dynamic_cast<PK_BDF_Default>(sub_pks_[0]);
  sed_transport_pk_ = sub_pks_[1];


  const Epetra_MultiVector& dz =
    *S->Get<CompositeVector>(elevation_increase_key_).ViewComponent("cell", false);

  dz_accumul_ = Teuchos::rcp(new Epetra_MultiVector(dz));
  dz_accumul_->PutScalar(0.);
}

void
Morphology_PK::CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S)
{
  sed_transport_pk_->CommitStep(t_old, t_new, S);

  S->set_time(t_new);

  Key elev_key = Keys::readKey(*plist_, domain_, "elevation", "elevation");
  Key slope_key = Keys::readKey(*plist_, domain_, "slope magnitude", "slope_magnitude");

  bool chg = S_->GetEvaluator(elev_key)->HasFieldChanged(S_.ptr(), elev_key);
  if (chg) {
    Teuchos::RCP<CompositeVector> elev = S->GetPtrW<CompositeVector>(elev_key, elev_key);
    Teuchos::RCP<CompositeVector> slope = S->GetPtrW<CompositeVector>(slope_key, slope_key);
    Teuchos::RCP<CompositeVector> vc = S->GetPtrW<CompositeVector>(vertex_coord_key_ss_, "state");
    Teuchos::RCP<CompositeVector> dz =
      S->GetPtrW<CompositeVector>(elevation_increase_key_, "state");

    *elev = *S_->GetPtrW<CompositeVector>(elev_key, elev_key);
    *slope = *S_->GetPtrW<CompositeVector>(slope_key, slope_key);
    *vc = *S_->GetPtrW<CompositeVector>(vertex_coord_key_ss_, "state");
  }

  Key biomass_key = Keys::getKey(domain_, "biomass");
  chg = S_->GetEvaluator(biomass_key)->HasFieldChanged(S_.ptr(), biomass_key);
  // if (chg)
}

// -----------------------------------------------------------------------------
// Advance each sub-PK individually, returning a failure as soon as possible.
// -----------------------------------------------------------------------------
bool
Morphology_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  bool fail = false;
  Teuchos::OSTab tab = vo_->getOSTab();

  double dt_step;
  dt_step = S_inter_->final_time() - S_inter_->initial_time();

  if (dt_step < dt_MPC_) {
    std::stringstream messagestream;
    messagestream << "Actual step is less than prescribed MPC time step.";
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }


  Teuchos::RCP<Field_Scalar> msl_rcp =
    Teuchos::rcp_dynamic_cast<Field_Scalar>(S_inter_->GetField("msl", "state"));
  msl_rcp->Compute(t_old);

  Key elev_key = Keys::readKey(*plist_, domain_, "elevation", "elevation");
  Epetra_MultiVector& dz =
    *S_next_->GetW<CompositeVector>(elevation_increase_key_, "state").ViewComponent("cell", false);
  dz.PutScalar(0.);

  flow_pk_->ResetTimeStepper(t_old);

  S_inter_->set_intermediate_time(t_old);
  S_next_->set_intermediate_time(t_old);
  double dt_done = 0;
  double dt_next = flow_pk_->get_dt();
  double t_DNS_end = t_old + dt_step / MSF_; // end of direct numerical simulation

  bool done = false;
  int ncycles = 0;

  while (!done) {
    dt_next = flow_pk_->get_dt();
    if (t_old + dt_done + dt_next > t_DNS_end) { dt_next = t_DNS_end - t_old - dt_done; }

    fail = true;
    while (fail) {
      S_next_->set_time(t_old + dt_done + dt_next);
      S_inter_->set_time(t_old + dt_done);
      fail = flow_pk_->AdvanceStep(t_old + dt_done, t_old + dt_done + dt_next, reinit);
      fail |= !flow_pk_->ValidStep();

      if (fail) {
        if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) *vo_->os() << "Master step is failed\n";
        dt_next = flow_pk_->get_dt();
      }
    }

    master_dt_ = dt_next;

    flow_pk_->CalculateDiagnostics(S_next_);
    flow_pk_->CommitStep(t_old + dt_done, t_old + dt_done + dt_next, S_next_);

    //S_next_->WriteStatistics(vo_);
    slave_dt_ = sed_transport_pk_->get_dt();
    if (slave_dt_ > master_dt_) slave_dt_ = master_dt_;
    if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
      *vo_->os() << "Slave dt=" << slave_dt_ << " Master dt=" << master_dt_ << "\n";

    fail = sed_transport_pk_->AdvanceStep(t_old + dt_done, t_old + dt_done + dt_next, reinit);

    if (fail) {
      dt_next /= 2;
    } else {
      S_inter_->set_intermediate_time(t_old + dt_done + dt_next);
      sed_transport_pk_->CommitStep(t_old + dt_done, t_old + dt_done + dt_next, S_next_);
      dt_done += dt_next;

      // we're done with this time step, copy the state
      *S_inter_ = *S_next_;
    }
    ncycles++;


    // check for done condition
    done = (std::abs(t_old + dt_done - t_DNS_end) / (t_DNS_end - t_old) <
            0.1 * min_dt_) ||   // finished the step
           (dt_next < min_dt_); // failed
  }

  double max_dz, min_dz;
  dz.MinValue(&min_dz);
  dz.MaxValue(&max_dz);
  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
    *vo_->os() << "min " << min_dz << " max " << max_dz << "\n";

  dz.Scale(MSF_);

  dz_accumul_->Update(1, dz, 1);
  Update_MeshVertices_(S_next_.ptr());


  bool chg = S_next_->GetEvaluator(elev_key)->HasFieldChanged(S_next_.ptr(), elev_key);
  Epetra_MultiVector& elev_cell =
    *S_next_->GetW<CompositeVector>(elev_key, elev_key).ViewComponent("cell", false);
  //S_next_ -> GetEvaluator("surface-slope)->HasFieldChanged(S_next_.ptr(), elev_key);


  return fail;


  // // advance the slave, subcycling if needed
  // S_next_->set_intermediate_time(t_old);
  // bool done = false;
  // double dt_next = slave_dt_;
  // double dt_done = 0.;
  // int ncycles = 0;

  // while (!done) {
  //   // do not overstep
  //   if (t_old + dt_done + dt_next > t_new) {
  //     dt_next = t_new - t_old - dt_done;
  //   }
  //   // take the step
  //   fail = sed_transport_pk_->AdvanceStep(t_old + dt_done, t_old + dt_done + dt_next, reinit);
  //   ncycles ++;
  //   if (fail) {
  //     // if fail, cut the step and try again
  //     dt_next /= 2;
  //   } else {
  //     // if success, commit the state and increment to next intermediate
  //     // -- etc: unclear if state should be commited or not?
  //     // set the intermediate time
  //     S_next_ -> set_intermediate_time(t_old + dt_done + dt_next);
  //     //S_next_ -> set_intermediate_time(t_old + dt_done + dt_next);
  //     sed_transport_pk_->CommitStep(t_old + dt_done, t_old + dt_done + dt_next, S_next_);
  //     //sed_transport_pk_->CommitStep(t_old + dt_done, t_old + dt_done + dt_next, S_next_);
  //     dt_done += dt_next;
  //   }
  //   // check for done condition
  //   done = (std::abs(t_old + dt_done - t_new) / (t_new - t_old) < 0.1*min_dt_) || // finished the step
  //       (dt_next  < min_dt_); // failed
  // }


  // if (std::abs(t_old + dt_done - t_new) / (t_new - t_old) < 0.1*min_dt_) {
  //   // done, success
  //   if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) *vo_->os()<<"Slave step is successful after "
  //                                                            <<ncycles <<" subcycles\n";
  //   fail = false;
  // } else {
  //   if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) *vo_->os()<<"Slave step is failed after "
  //                                                            <<ncycles <<" subcycles\n";
  //   fail = true;
  // }

  return fail;
}

void
Morphology_PK::Initialize_MeshVertices_(const Teuchos::Ptr<State>& S,
                                        Teuchos::RCP<const AmanziMesh::Mesh> mesh,
                                        Key vert_field_key)
{
  // spatial dimension
  int dim = mesh->getSpaceDimension();
  Amanzi::AmanziGeometry::Point coords(dim);
  // number of vertices
  int nV = mesh->getNumEntities(Amanzi::AmanziMesh::Entity_kind::NODE, Amanzi::AmanziMesh::Parallel_kind::OWNED);

  Epetra_MultiVector& vc =
    *S->GetPtrW<CompositeVector>(vert_field_key, "state")->ViewComponent("node", false);

  // search the id of the mid point on the top
  for (int iV = 0; iV < nV; iV++) {
    // get the coords of the node
    coords = mesh->getNodeCoordinate(iV);
    for (int s = 0; s < dim; ++s) { vc[s][iV] = coords[s]; }
  }

  S->GetPtrW<CompositeVector>(vert_field_key, "state")->ScatterMasterToGhosted("node");
  S->GetField(vert_field_key, "state")->set_initialized();
}

void
Morphology_PK::Update_MeshVertices_(const Teuchos::Ptr<State>& S)
{
  // spatial dimension
  int dim = mesh_ss_->getSpaceDimension();
  Amanzi::AmanziGeometry::Point coords(dim);
  // number of vertices

  Epetra_MultiVector& vc =
    *S->GetPtrW<CompositeVector>(vertex_coord_key_ss_, "state")->ViewComponent("node", true);


  const Epetra_MultiVector& dz =
    *S->Get<CompositeVector>(elevation_increase_key_).ViewComponent("cell");

  int ncells = dz.MyLength();

  AmanziMesh::Entity_ID_List nodes, cells;
  double xyz[3];

  for (int c = 0; c < ncells; c++) {
    AmanziMesh::Entity_ID domain_face;
    domain_face = mesh_->getEntityParent(AmanziMesh::Entity_kind::CELL, c);

    nodes = mesh_ss_->getFaceNodes(domain_face);
    int nnodes = nodes.size();
    for (int i = 0; i < nnodes; i++) {
      coords = mesh_ss_->getNodeCoordinate(nodes[i]);
      cells = mesh_ss_->getNodeCells(nodes[i], Amanzi::AmanziMesh::Parallel_kind::OWNED);
      int nsize = cells.size();
      double old = coords[2];

      coords[2] += dz[0][c] / nsize;
      vc[2][nodes[i]] += dz[0][c] / nsize;

      // coords[2] += 0.1 / nsize;
      // vc[2][nodes[i]] += 0.1 / nsize;

      // if (c==0){
      //   std::cout<<"COORDS "<<i<<" "<<old <<" "<<coords[2]<<" "<<vc[2][nodes[i]]<<"\n";
      // }

      mesh_ss_->node_set_coordinates(nodes[i], coords);
    }
  }

  deform_eval_ =
    Teuchos::rcp_dynamic_cast<EvaluatorPrimary>(S->GetEvaluator(elevation_increase_key_));
  deform_eval_->SetChanged(S.ptr());


  S->GetPtrW<CompositeVector>(vertex_coord_key_, "state")->ScatterMasterToGhosted("node");
}


} // namespace Amanzi
