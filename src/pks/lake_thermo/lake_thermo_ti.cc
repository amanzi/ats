/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon
------------------------------------------------------------------------- */

#include "boost/math/special_functions/fpclassify.hpp"
#include "boost/test/floating_point_comparison.hpp"

#include "Debugger.hh"
#include "BoundaryFunction.hh"
#include "FieldEvaluator.hh"
#include "lake_thermo_pk.hh"
#include "Op.hh"

namespace Amanzi {
namespace LakeThermo {

using namespace ATS::Mesh;

#define DEBUG_FLAG 1
#define MORE_DEBUG_FLAG 0

// Lake_Thermo_PK is a BDFFnBase
// -----------------------------------------------------------------------------
// computes the non-linear functional g = g(t,u,udot)
// -----------------------------------------------------------------------------
void Lake_Thermo_PK::FunctionalResidual(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
    Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g) {
  Teuchos::OSTab tab = vo_->getOSTab();

  bool ice_cover_ = false; // first always assume that there is no ice

  // get temperature
  Teuchos::RCP<const CompositeVector> temp = S_inter_->GetFieldData(temperature_key_);

  Epetra_MultiVector& g_c = *g->Data()->ViewComponent("cell",false);
  unsigned int ncomp = g_c.MyLength();

  const Epetra_MultiVector& temp_v = *temp->ViewComponent("cell",false);

  // cell volumes
  const Epetra_MultiVector& cv =
         *S_inter_->GetFieldData(Keys::getKey(domain_,"cell_volume"))->ViewComponent("cell",false);

  // get conductivity
  S_inter_->GetFieldEvaluator(conductivity_key_)->HasFieldChanged(S_inter_.ptr(), name_);
  const Epetra_MultiVector& lambda_c =
      *S_inter_->GetFieldData(conductivity_key_)->ViewComponent("cell",false);

  // increment, get timestep
  niter_++;
  double h = t_new - t_old;

  // pointer-copy temperature into states and update any auxilary data
  Solution_to_State(*u_new, S_next_);
  Teuchos::RCP<CompositeVector> u = u_new->Data();

#if DEBUG_FLAG
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "----------------------------------------------------------------" << std::endl
    << "Residual calculation: t0 = " << t_old
    << " t1 = " << t_new << " h = " << h << std::endl;

  // dump u_old, u_new
  db_->WriteCellInfo(true);
  std::vector<std::string> vnames;
  vnames.push_back("T_old"); vnames.push_back("T_new");
  std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
  vecs.push_back(S_inter_->GetFieldData(key_).ptr()); vecs.push_back(u.ptr());
  db_->WriteVectors(vnames, vecs, true);

  // vnames[0] = "sl"; vnames[1] = "si";
  // vecs[0] = S_next_->GetFieldData("saturation_liquid").ptr();
  // vecs[1] = S_next_->GetFieldData("saturation_ice").ptr();
  // db_->WriteVectors(vnames, vecs, false);
#endif

  // update boundary conditions
  bc_temperature_->Compute(t_new);
  bc_diff_flux_->Compute(t_new);
  //  bc_flux_->Compute(t_new);
  UpdateBoundaryConditions_(S_next_.ptr());

  S_inter_->GetFieldData(depth_key_,name_)->PutScalar(h_);

  /*
  // generate new mesh

  Comm_ptr_type comm = Amanzi::getDefaultComm();

  // create a mesh framework
//  Teuchos::ParameterList regions_list = plist_->get<Teuchos::ParameterList>("regions");
//  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, regions_list, *comm));
  Teuchos::RCP<AmanziGeometry::GeometricModel> gm = Teuchos::rcp(new AmanziGeometry::GeometricModel(3, plist_->sublist("regions"), *comm));

  Teuchos::ParameterList& meshes_list = plist_->sublist("mesh");
  Teuchos::ParameterList gen_mesh_plist;
  gen_mesh_plist = meshes_list.sublist("lake scaled");
  double high_coords[3] = {1.,1.,h_};
  gen_mesh_plist.set("domain high coordinate", high_coords);

//  mesh_->set_parameter_list(gen_mesh_plist);

  // create a mesh
  bool request_faces = true, request_edges = false;
  Amanzi::AmanziMesh::MeshFactory meshfactory(comm,gm);
  meshfactory.set_preference(Amanzi::AmanziMesh::Preference({Amanzi::AmanziMesh::Framework::MSTK}));

  mesh_scaled_ = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, h_, 1, 1, ncomp, request_faces, request_edges);
  */

  /*
  Comm_ptr_type comm;
  Teuchos::RCP<AmanziGeometry::GeometricModel> gm;
  Amanzi::State S;

  comm = getDefaultComm();
  gm = Teuchos::rcp(new AmanziGeometry::GeometricModel(3, plist_->sublist("regions"), *comm));
//  S = Teuchos::rcp(new State(plist_->sublist("state")));

  Teuchos::ParameterList& meshes_list = plist_->sublist("mesh");
  Teuchos::ParameterList gen_mesh_plist;
  gen_mesh_plist = meshes_list.sublist("lake scaled");
  double high_coords[3] = {1.,1.,h_};
  gen_mesh_plist.set("domain high coordinate", high_coords);

  VerboseObject vo(comm, "ATS Mesh Factory", meshes_list);
//  ATS::Mesh::createMesh(meshes_list.sublist("lake scaled"), comm, gm, S, vo);
 * */

  // zero out residual
  Teuchos::RCP<CompositeVector> res = g->Data();
  res->PutScalar(0.0);

  // diffusion term, implicit
  ApplyDiffusion_(S_next_.ptr(), res.ptr());
#if DEBUG_FLAG
  db_->WriteVector("K",S_next_->GetFieldData(conductivity_key_).ptr(),true);
  db_->WriteVector("res (diff)", res.ptr(), true);
#endif

  // source terms
  AddSources_(S_next_.ptr(), res.ptr());
#if DEBUG_FLAG
  db_->WriteVector("res (src)", res.ptr());
#endif

  // accumulation term
  AddAccumulation_(res.ptr());
#if DEBUG_FLAG
  vnames[0] = "T_old";
  vnames[1] = "T_new";
  vecs[0] = S_inter_->GetFieldData(temperature_key_).ptr();
  vecs[1] = S_next_->GetFieldData(temperature_key_).ptr();
  db_->WriteVectors(vnames, vecs, true);
  db_->WriteVector("res (acc)", res.ptr());
#endif

  // advection term
  if (implicit_advection_) {
    AddAdvection_(S_next_.ptr(), res.ptr(), true);
  } else {
    AddAdvection_(S_inter_.ptr(), res.ptr(), true);
  }
#if DEBUG_FLAG
  db_->WriteVector("res (adv)", res.ptr());
#endif

  // Dump residual to state for visual debugging.
#if MORE_DEBUG_FLAG
  if (niter_ < 23) {
    std::stringstream namestream;
    namestream << domain_prefix_ << "energy_residual_" << niter_;
    *S_next_->GetFieldData(namestream.str(),name_) = *res;

    std::stringstream solnstream;
    solnstream << domain_prefix_ << "energy_solution_" << niter_;
    *S_next_->GetFieldData(solnstream.str(),name_) = *u;
  }
#endif

// const Epetra_MultiVector& g_c = *res->ViewComponent("cell", false);

// unsigned int ncells = g_c.MyLength();
// for (unsigned int c=0; c!=ncomp; ++c) {
//   g_c[0][c] = g_c[0][c]*cv[0][0];
// }

};


// -----------------------------------------------------------------------------
// Apply the preconditioner to u and return the result in Pu.
// -----------------------------------------------------------------------------
int Lake_Thermo_PK::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) {
#if DEBUG_FLAG
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "Precon application:" << std::endl;
  db_->WriteVector("T_res", u->Data().ptr(), true);
#endif

  // apply the preconditioner
  int ierr = preconditioner_->ApplyInverse(*u->Data(), *Pu->Data());

#if DEBUG_FLAG
  db_->WriteVector("PC*T_res", Pu->Data().ptr(), true);
#endif

  Pu->Data()->ViewComponent("boundary_face")->PutScalar(0.0); // correction 01/22/21, for correct boundary conditions

  return (ierr > 0) ? 0 : 1;
};

// -----------------------------------------------------------------------------
// Update the preconditioner at time t and u = up
// -----------------------------------------------------------------------------
void Lake_Thermo_PK::UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h) {
  // VerboseObject stuff.
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "Precon update at t = " << t << std::endl;

  // update state with the solution up.

  AMANZI_ASSERT(std::abs(S_next_->time() - t) <= 1.e-4*t);
  PK_PhysicalBDF_Default::Solution_to_State(*up, S_next_);

  Teuchos::RCP<const CompositeVector> temp = S_next_ -> GetFieldData(key_);

  // update boundary conditions
  bc_temperature_->Compute(S_next_->time());
  bc_diff_flux_->Compute(S_next_->time());
  //  bc_flux_->Compute(S_next_->time());
  UpdateBoundaryConditions_(S_next_.ptr());

  // div K_e grad u
  UpdateConductivityData_(S_next_.ptr());
  if (jacobian_) UpdateConductivityDerivativeData_(S_next_.ptr());

  Teuchos::RCP<const CompositeVector> conductivity =
      S_next_->GetFieldData(uw_conductivity_key_);

  // jacobian term
  Teuchos::RCP<const CompositeVector> dKdT = Teuchos::null;
  if (jacobian_) {
    if (!duw_conductivity_key_.empty()) {
      dKdT = S_next_->GetFieldData(duw_conductivity_key_);
    } else {
      dKdT = S_next_->GetFieldData(dconductivity_key_);
    }
  }

  // create local matrices
  preconditioner_->Init();
  preconditioner_diff_->SetScalarCoefficient(conductivity, dKdT);
  preconditioner_diff_->UpdateMatrices(Teuchos::null, temp.ptr());
  preconditioner_diff_->ApplyBCs(true, true, true);

  if (jacobian_) {
    Teuchos::RCP<CompositeVector> flux = Teuchos::null;

    flux = S_next_->GetFieldData(energy_flux_key_, name_);
    preconditioner_diff_->UpdateFlux(up->Data().ptr(), flux.ptr());

    preconditioner_diff_->UpdateMatricesNewtonCorrection(flux.ptr(), up->Data().ptr());
  }

  // update with accumulation terms
  // -- update the accumulation derivatives, de/dT
  S_next_->GetFieldEvaluator(energy_key_)
          ->HasFieldDerivativeChanged(S_next_.ptr(), name_, key_);
  const Epetra_MultiVector& de_dT = *S_next_->GetFieldData(Keys::getDerivKey(energy_key_, key_))
          ->ViewComponent("cell",false);
  unsigned int ncells = de_dT.MyLength();

  CompositeVector acc(S_next_->GetFieldData(energy_key_)->Map());
  auto& acc_c = *acc.ViewComponent("cell", false);

#if DEBUG_FLAG
  //  db_->WriteVector("    de_dT", S_next_->GetFieldData(Keys::getDerivKey(temperature_key_, key_)).ptr());
#endif

  if (coupled_to_subsurface_via_temp_ || coupled_to_subsurface_via_flux_) {
    // do not add in de/dT if the height is 0
    const Epetra_MultiVector& pres = *S_next_->GetFieldData(Keys::getKey(domain_,"pressure"))
            ->ViewComponent("cell",false);
    const double& patm = *S_next_->GetScalarData("atmospheric_pressure");

    for (unsigned int c=0; c!=ncells; ++c) {
      acc_c[0][c] = pres[0][c] >= patm ? de_dT[0][c] / h : 0.;
    }
  } else {
    if (precon_used_) {
      for (unsigned int c=0; c!=ncells; ++c) {
        //      AMANZI_ASSERT(de_dT[0][c] > 1.e-10);
        // ?? Not using e_bar anymore apparently, though I didn't think we were ever.  Need a nonzero here to ensure not singlar.
        acc_c[0][c] = std::max(de_dT[0][c], 1.e-12) / h;
      }
    } else {
      for (unsigned int c=0; c!=ncells; ++c) {
        //      AMANZI_ASSERT(de_dT[0][c] > 1.e-10);
        // ?? Not using e_bar anymore apparently, though I didn't think we were ever.  Need a nonzero here to ensure not singlar.
        // apply a diagonal shift manually for coupled problems
        acc_c[0][c] = de_dT[0][c] / h + 1.e-6;
      }
    }
  }
  preconditioner_acc_->AddAccumulationTerm(acc, "cell");

  // -- update preconditioner with source term derivatives if needed
  AddSourcesToPrecon_(S_next_.ptr(), h);

  //  // update with advection terms
  //  if (implicit_advection_ && implicit_advection_in_pc_) {
  //    Teuchos::RCP<const CompositeVector> mass_flux = S_next_->GetFieldData(flux_key_);
  //    S_next_->GetFieldEvaluator(enthalpy_key_)
  //        ->HasFieldDerivativeChanged(S_next_.ptr(), name_, key_);
  ////    Teuchos::RCP<const CompositeVector> dhdT = S_next_->GetFieldData(Keys::getDerivKey(enthalpy_key_, key_));
  //    Teuchos::RCP<const CompositeVector> dhdT = S_next_->GetFieldData(Keys::getDerivKey(temperature_key_, key_));
  //    preconditioner_adv_->Setup(*mass_flux);
  //    preconditioner_adv_->SetBCs(bc_adv_, bc_adv_);
  //    preconditioner_adv_->UpdateMatrices(mass_flux.ptr(), dhdT.ptr());
  //    ApplyDirichletBCsToEnthalpy_(S_next_.ptr());  //!!!!!!!!!!!!! IMPORTANT !!!!!!
  //    preconditioner_adv_->ApplyBCs(false, true, false);
  //
  //  }

  // Apply boundary conditions.
  preconditioner_diff_->ApplyBCs(true, true, true);
};

double Lake_Thermo_PK::ErrorNorm(Teuchos::RCP<const TreeVector> u,
    Teuchos::RCP<const TreeVector> du)
{
  Teuchos::OSTab tab = vo_->getOSTab();

  // Relative error in cell-centered temperature
  const Epetra_MultiVector& uc = *u->Data()->ViewComponent("cell", false);
  const Epetra_MultiVector& duc = *du->Data()->ViewComponent("cell", false);

  int ncells_owned = mesh_->num_entities(Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED);

  double error_t(0.0);
  double ref_temp(273.0);
  for (int c = 0; c < ncells_owned; c++) {
    double tmp = fabs(duc[0][c]) / (fabs(uc[0][c] - ref_temp) + ref_temp);
    if (tmp > error_t) {
      error_t = tmp;
    }
  }

  // Cell error is based upon error in energy conservation relative to
  // a characteristic energy
  double error_e(0.0);

  double error = std::max(error_t, error_e);

#ifdef HAVE_MPI
  double buf = error;
  du->Data()->Comm()->MaxAll(&buf, &error, 1);  // find the global maximum
#endif

  return error;
}

} // namespace LakeThermo
} // namespace Amanzi
