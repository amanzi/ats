
#include <iostream>

#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include "Epetra_SerialComm.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

#include "VerboseObject_objs.hh"
#include "VerboseObject.hh"

#include "ats_version.hh"
#include "amanzi_version.hh"
#include "tpl_versions.h"

#include "AmanziComm.hh"
#include "AmanziTypes.hh"
#include "GeometricModel.hh"
#include "State.hh"

#include "errors.hh"
#include "exceptions.hh"
#include "dbc.hh"

#include "ats_mesh_factory.hh"
#include "elm_ats_coordinator.hh"
#include "elm_ats_driver.hh"

// registration files
#include "state_evaluators_registration.hh"

#include "ats_relations_registration.hh"
#include "ats_transport_registration.hh"
#include "ats_energy_pks_registration.hh"
#include "ats_energy_relations_registration.hh"
#include "ats_flow_pks_registration.hh"
#include "ats_flow_relations_registration.hh"
#include "ats_deformation_registration.hh"
#include "ats_bgc_registration.hh"
#include "ats_surface_balance_registration.hh"
#include "ats_mpc_registration.hh"
//#include "ats_sediment_transport_registration.hh"
#include "mdm_transport_registration.hh"
#include "multiscale_transport_registration.hh"
#ifdef ALQUIMIA_ENABLED
#include "pks_chemistry_registration.hh"
#endif

// include fenv if it exists
#include "boost/version.hpp"
#if (BOOST_VERSION / 100 % 1000 >= 46)
#include "boost/config.hpp"
#ifndef BOOST_NO_FENV_H
#ifdef _GNU_SOURCE
#define AMANZI_USE_FENV
#include "boost/detail/fenv.hpp"
#endif
#endif
#endif

#include "boost/filesystem.hpp"

#include "pk_helpers.hh"

namespace ATS {

void
ELM_ATSDriver::setup(MPI_Fint *f_comm, const char *infile)
{
  // -- create communicator & get process rank
  //auto comm = Amanzi::getDefaultComm();
  auto c_comm = MPI_Comm_f2c(*f_comm);
  auto comm = setComm(c_comm);
  auto rank = comm->MyPID();

  // convert input file to std::string for easier handling
  // infile must be null-terminated
  std::string input_filename(infile);

  // parse the input file and check validity
  if (input_filename.empty()) {
    if (rank == 0)
      std::cerr << "ERROR: no input file provided" << std::endl;
  } else if (!boost::filesystem::exists(input_filename)) {
    if (rank == 0)
      std::cerr << "ERROR: input file \"" << input_filename << "\" does not exist." << std::endl;
  }

  // -- parse input file
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(input_filename);

  // -- set default verbosity level to no output
  Amanzi::VerboseObject::global_default_level = Teuchos::VERB_NONE;

  // create the geometric model and regions
  Teuchos::ParameterList reg_params = plist->sublist("regions");
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
    Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, reg_params, *comm) );

  // Create the state.
  Teuchos::ParameterList state_plist = plist->sublist("state");
  S_ = Teuchos::rcp(new Amanzi::State(state_plist));

  // create and register meshes
  ATS::Mesh::createMeshes(*plist, comm, gm, *S_);

  // keys
  domain_sub_ = plist->get<std::string>("domain name", "domain");
  domain_srf_ = Amanzi::Keys::readDomainHint(*plist, domain_sub_, "subsurface", "surface");
  sub_src_key_ = Amanzi::Keys::readKey(*plist, domain_sub_, "subsurface source", "source_sink");
  srf_src_key_ = Amanzi::Keys::readKey(*plist, domain_srf_, "surface source", "source_sink");
  pres_key_ = Amanzi::Keys::readKey(*plist, domain_sub_, "pressure", "pressure");
  satl_key_ = Amanzi::Keys::readKey(*plist, domain_sub_, "saturation_liquid", "saturation_liquid");
  por_key_ = Amanzi::Keys::readKey(*plist, domain_sub_, "porosity", "porosity");


  srf_mol_dens_key_ = Amanzi::Keys::readKey(*plist, domain_srf_, "molar density", "molar_density_liquid");
  srf_mass_dens_key_ = Amanzi::Keys::readKey(*plist, domain_srf_, "mass density", "mass_density_liquid");


  //sub_mol_dens_key_
  //sub_mass_dens_key_

  // assume for now that mesh info has been communicated
  mesh_subsurf_ = S_->GetMesh(domain_sub_);
  mesh_surf_ = S_->GetMesh(domain_srf_);

  // build columns to allow indexing by column
  mesh_subsurf_->build_columns();

  // check that number of surface cells = number of columns
  ncolumns_ = mesh_surf_->num_entities(Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED);
  AMANZI_ASSERT(ncolumns_ == mesh_subsurf_->num_columns(false));

  // get num cells per column - include consistency check later
  auto& col_zero = mesh_subsurf_->cells_of_column(0);
  std::size_t ncol_cells_ = col_zero.size();

  // require primary variables
  // -- subsurface water source
  S_->Require<Amanzi::CompositeVector,Amanzi::CompositeVectorSpace>(sub_src_key_, Amanzi::Tags::NEXT,  sub_src_key_)
    .SetMesh(mesh_subsurf_)->SetComponent("cell", Amanzi::AmanziMesh::CELL, 1);
  RequireEvaluatorPrimary(sub_src_key_, Amanzi::Tags::NEXT, *S_);
  // -- surface water source-sink
  S_->Require<Amanzi::CompositeVector,Amanzi::CompositeVectorSpace>(srf_src_key_, Amanzi::Tags::NEXT,  srf_src_key_)
    .SetMesh(mesh_surf_)->SetComponent("cell", Amanzi::AmanziMesh::CELL, 1);
  RequireEvaluatorPrimary(srf_src_key_, Amanzi::Tags::NEXT, *S_);
  // -- porosity
  //S_->Require<Amanzi::CompositeVector,Amanzi::CompositeVectorSpace>(por_key_, Amanzi::Tags::NEXT,  por_key_)
  //  .SetMesh(mesh_subsurf_)->SetComponent("cell", Amanzi::AmanziMesh::CELL, 1);
  //RequireEvaluatorPrimary(por_key_, Amanzi::Tags::NEXT, *S_);

  // create ELM coordinator object
  elm_coordinator_ = std::make_unique<ELM_ATSCoordinator>(*plist, S_, comm);
  // call coordinator setup
  elm_coordinator_->setup();
}

void ELM_ATSDriver::initialize()
{
  elm_coordinator_->initialize();

  S_->GetW<Amanzi::CompositeVector>(sub_src_key_, Amanzi::Tags::NEXT, sub_src_key_).PutScalar(0.);
  S_->GetRecordW(sub_src_key_, Amanzi::Tags::NEXT, sub_src_key_).set_initialized();
  S_->GetW<Amanzi::CompositeVector>(srf_src_key_, Amanzi::Tags::NEXT, srf_src_key_).PutScalar(0.);
  S_->GetRecordW(srf_src_key_, Amanzi::Tags::NEXT, srf_src_key_).set_initialized();
  //S_->GetW<Amanzi::CompositeVector>(por_key_, Amanzi::Tags::NEXT, por_key_).PutScalar(0.);
  //S_->GetRecordW(por_key_, Amanzi::Tags::NEXT, por_key_).set_initialized();

}


void ELM_ATSDriver::advance(double *dt)
{
  // Get incoming state from ELM
  Epetra_MultiVector& srf_water_src = *S_->GetW<Amanzi::CompositeVector>(srf_src_key_, Amanzi::Tags::NEXT, srf_src_key_)
    .ViewComponent("cell", false);
  Epetra_MultiVector& sub_water_src = *S_->GetW<Amanzi::CompositeVector>(sub_src_key_, Amanzi::Tags::NEXT, sub_src_key_)
    .ViewComponent("cell", false);
  //Epetra_MultiVector& porosity = *S_->GetW<Amanzi::CompositeVector>(por_key_, Amanzi::Tags::NEXT, por_key_)
  //  .ViewComponent("cell", false);
  ChangedEvaluatorPrimary(srf_src_key_, Amanzi::Tags::NEXT, *S_);
  ChangedEvaluatorPrimary(sub_src_key_, Amanzi::Tags::NEXT, *S_);
  //ChangedEvaluatorPrimary(por_key_, Amanzi::Tags::NEXT, *S_);




  const Epetra_MultiVector& testmolar = *S_->Get<Amanzi::CompositeVector>("molar_density_liquid", Amanzi::Tags::NEXT)
    .ViewComponent("cell", false);
  const Epetra_MultiVector& testmass = *S_->Get<Amanzi::CompositeVector>("mass_density_liquid", Amanzi::Tags::NEXT)
    .ViewComponent("cell", false);
//
for (Amanzi::AmanziMesh::Entity_ID col=0; col!=ncolumns_; ++col) {
  auto& col_iter = mesh_subsurf_->cells_of_column(col);
  for (std::size_t i=0; i!=col_iter.size(); ++i) {
    std::cout << "DENSITIES:::  " << testmolar[0][col_iter[i]] << "   " << testmass[0][col_iter[i]] << std::endl;

  }
}





  elm_coordinator_->advance(*dt);

  // update ATS->ELM data if necessary
  S_->GetEvaluator(pres_key_, Amanzi::Tags::NEXT).Update(*S_, pres_key_);
  const Epetra_MultiVector& pres = *S_->Get<Amanzi::CompositeVector>(pres_key_, Amanzi::Tags::NEXT).ViewComponent("cell", false);
  S_->GetEvaluator(satl_key_, Amanzi::Tags::NEXT).Update(*S_, satl_key_);
  const Epetra_MultiVector& satl = *S_->Get<Amanzi::CompositeVector>(satl_key_, Amanzi::Tags::NEXT).ViewComponent("cell", false);
  //S_->GetEvaluator(por_key_, Amanzi::Tags::NEXT).Update(*S_, por_key_);
  //const Epetra_MultiVector& poro = *S_->Get<Amanzi::CompositeVector>(por_key_, Amanzi::Tags::NEXT).ViewComponent("cell", false);

}

void ELM_ATSDriver::advance_test()
{
  // use dt from ATS for now
  double dt = elm_coordinator_->get_dt(false);

  while (S_->get_time() < elm_coordinator_->get_end_time()) {
    advance(&dt);
    dt = elm_coordinator_->get_dt(false);
  };
}

/* assume that incoming data is in form
soil_infiltration(ncolumns)
soil_evaporation(ncolumns)
root_transpiration(ncells)



*/

void ELM_ATSDriver::set_sources(double* soil_infiltration, double* soil_evaporation, int *ncols)
{
  Epetra_MultiVector& surf_ss = *S_->GetW<Amanzi::CompositeVector>(srf_src_key_, Amanzi::Tags::NEXT, srf_src_key_)
      .ViewComponent("cell", false);
  const Epetra_MultiVector& mol_dens = *S_->Get<Amanzi::CompositeVector>(srf_mol_dens_key_, Amanzi::Tags::NEXT)
      .ViewComponent("cell", false);
  const Epetra_MultiVector& mass_dens = *S_->Get<Amanzi::CompositeVector>(srf_mass_dens_key_, Amanzi::Tags::NEXT)
      .ViewComponent("cell", false);
  
  AMANZI_ASSERT(*ncols == ncolumns_ == surf_ss.MyLength());
  //AMANZI_ASSERT(*ncells == subsurf_ss.MyLength());

  for (Amanzi::AmanziMesh::Entity_ID col=0; col!=ncolumns_; ++col) {
    double mol_h20_kg = mol_dens[0][col] / mass_dens[0][col];
    surf_ss[0][col] = soil_evaporation[col] * mol_h20_kg;
    std::cout << "surf_ss1::   " << surf_ss[0][col] << "  " soil_evaporation[col] << "  " <<  mol_h20_kg << std::endl;
    surf_ss[0][col] += (soil_infiltration[col] * mol_h20_kg);
    std::cout << "surf_ss2::   " << surf_ss[0][col] << "  " soil_infiltration[col] << "  " <<  mol_h20_kg << std::endl;

  }

}


// helper function for pushing ELM data to ats columns
// assumes incoming fortran array was in form arr(ncells, ncolumns)
// where ncells dimension varies fastest 
//void ELM_ATSDriver::push_to_columns(Epetra_Vector& ats_vec, double* elm_data)
//{
//  for (Amanzi::AmanziMesh::Entity_ID col=0; col!=ncolumns_; ++col) {
//    auto& col_iter = mesh_subsurf_->cells_of_column(col);
//    for (std::size_t i=0; i!=col_iter.size(); ++i) {
//      ats_vec[col_iter[i]] = elm_data[col*ncol_cells_+i];
//    }
//  }
//}

// helper function for pushing field to column
// taken from fates_pk
//void ELM_ATSDriver::push_to_column(Amanzi::AmanziMesh::Entity_ID col, const Epetra_Vector& vec,
//                               double* col_vec, int ncol) {
//
//  auto& col_iter = mesh_subsurf_->cells_of_column(col);
//  for (std::size_t i=0; i!=col_iter.size(); ++i) {
//    col_vec[i] = vec[col_iter[i]];
//  }
//}



} // namespace ATS
