/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
  License: see $ATS_DIR/COPYRIGHT
  Author: 
*/

/*!

Requires the following dependencies:


*/

#include "Key.hh"
#include "distributed_tiles_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {


DistributedTilesRateEvaluator::DistributedTilesRateEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist),
    compatibility_checked_(false)
{


  domain_ = Keys::getDomain(my_key_);

  subsurface_marks_key_ = Keys::readKey(plist, domain_, "catchments_id", "catchments_id");
  sources_key_ = plist.get<std::string>("accumulated source key", "subdomain_sources");
  factor_key_ = plist.get<std::string>("factor field key", "");
  num_component_ = plist.get<int>("number of components", 1);
    
  // surface_marks_key_ = Keys::readKey(plist, domain_surf_, "ditch marks", "ditch_marks");
  // surf_len_key_ = Keys::readKey(plist, domain_surf_, "ditch length", "ditch_length");
  
  num_ditches_ = plist.get<int>("number of ditches");

  implicit_ = plist.get<bool>("implicit drainage", true);

  dependencies_.insert(subsurface_marks_key_);
  pres_key_ = Keys::readKey(plist, domain_, "pressure", "pressure");
  mol_dens_key_ = Keys::readKey(plist, domain_, "molar density", "molar_density_liquid");
  mass_dens_key_ = Keys::readKey(plist, domain_, "mass density", "mass_density_liquid");
  visc_key_ = Keys::readKey(plist, domain_, "viscosity liquid", "viscosity_liquid");
  dependencies_.insert(pres_key_);
  dependencies_.insert(mol_dens_key_);
  dependencies_.insert(visc_key_);
  dependencies_.insert(mass_dens_key_);

  if (factor_key_ != "") dependencies_.insert(factor_key_);
  
  times_.resize(2);
  times_[0] = -1.0; times_[1] = -1.0;

  ka_ = plist.get<double>("permeability above [m^2]");
  kb_ = plist.get<double>("permeability below [m^2]");  
  d_ = plist.get<double>("equivalent distance to impermeable boundary [m]");
  L_ = plist.get<double>("drain spacing [m]");
  th_ = plist.get<double>("drain layer thickness [m]", -1.0);
  p_enter_ = plist.get<double>("entering pressure [Pa]", 101325);

}

// Required methods from SecondaryVariableFieldEvaluator
void
DistributedTilesRateEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
                 const Teuchos::Ptr<CompositeVector>& result)
{

  auto mesh = S->GetMesh(domain_);
  double t0 = S->time();
  double dt = *S->GetScalarData("dt");
  //  double t1 = S->final_time();


  Teuchos::RCP<Field> src_field =  S->GetField(sources_key_, "state");

  Teuchos::RCP<const Epetra_Vector> gvec = S->GetConstantVectorData("gravity");
  int z_index = mesh->space_dimension() - 1; 
  double gr = -(*gvec)[z_index];

  const auto& mass_dens = *S->GetFieldData(mass_dens_key_)->ViewComponent("cell", false);
  const auto& visc = *S->GetFieldData(visc_key_)->ViewComponent("cell", false);
  
  const auto& pres = *S->GetFieldData(pres_key_)->ViewComponent("cell", false);
  const auto& dens = *S->GetFieldData(mol_dens_key_)->ViewComponent("cell", false);
  const auto& cv =
      *S->GetFieldData(Keys::getKey(domain_,"cell_volume"))->ViewComponent("cell",false);
  const auto& sub_marks = *S->GetFieldData(subsurface_marks_key_)->ViewComponent("cell", false);
  
  auto& sub_sink = *result->ViewComponent("cell",false);
  sub_sink.PutScalar(0);
  
  AmanziMesh::Entity_ID ncells = sub_marks.MyLength();
  Teuchos::RCP<Epetra_Vector> src_vec = src_field->GetConstantVectorData();

  double total = 0.0;
  src_vec->PutScalar(0.0);
  sub_sink.PutScalar(0.0);
  int num_vectors = 1;
  int test = sub_sink.NumVectors();

  if (factor_key_!=""){
    num_vectors = S->GetFieldData(factor_key_)->ViewComponent("cell",false)->NumVectors();
    AMANZI_ASSERT(num_vectors == sub_sink.NumVectors());
    AMANZI_ASSERT(num_vectors == num_component_);
  }

  if (abs(dt) > 1e-13) {
    for (AmanziMesh::Entity_ID c=0; c!=ncells; ++c) {
      if (sub_marks[0][c] > 0) {
        if (th_ < 0){
          AmanziMesh::Entity_ID_List faces;
          std::vector<int> dirs;
          mesh->cell_get_faces_and_dirs(c, &faces, &dirs);
          double zmax = -1e+98;
          double zmin =  1e+98;          
          for (auto f : faces){
            auto xf=mesh->face_centroid(f);
            zmax = std::max(zmax, xf[z_index]);
            zmin = std::min(zmin, xf[z_index]);              
          }
          th_ = zmax - zmin;
          AMANZI_ASSERT(th_ > 1e-12);
        }

        
        double diff_p = std::min(p_enter_ - 1.5*pres[0][c], 0.0);
        double val = diff_p * dens[0][c] * (8.*kb_*d_)/(L_*L_*visc[0][c] *th_); //linear term [mol/(m^3 s)]
        val += diff_p * diff_p * dens[0][c] * (ka_/ (visc[0][c] * mass_dens[0][c] * L_ * L_ * gr * th_)); // nonlinear term [mol/(m^3 s)]

        if (factor_key_ != "") {
          const auto& factor = *S->GetFieldData(factor_key_)->ViewComponent("cell", false);
          for (int i = 0; i < num_component_; ++i) {
            sub_sink[i][c] = factor[i][c] * val;
            (*src_vec)[sub_marks[0][c] - 1 + i * num_ditches_] += factor[i][c] * val * dt * cv[0][c];
            total = total + factor[i][c] * val * dt * cv[0][c];
          }
        }else{
          sub_sink[0][c] = val;
          (*src_vec)[sub_marks[0][c] - 1] +=  val * dt * cv[0][c];
          total = total + val * dt * cv[0][c];
        }
      }
    }
  }
  // std::cout <<"Sink vector\n";
  // std::cout << *src_vec<<"\n";
  //std::cout<<"Total sink "<<my_key_<<" "<<total<<"\n";
  // total = 0.;
  // for (AmanziMesh::Entity_ID c=0; c!=ncells; ++c) {
  //   total += sub_sink[0][c] * cv[0][c];
  // }
  // std::cout<<"Total sink field "<<my_key_<<" "<<total<<"\n";

  
  Teuchos::RCP<const Comm_type> comm_p = S->GetMesh(domain_)->get_comm();
  Teuchos::RCP<const MpiComm_type> mpi_comm_p =
    Teuchos::rcp_dynamic_cast<const MpiComm_type>(comm_p);
  const MPI_Comm& comm = mpi_comm_p->Comm();

  double *src_vec_ptr;
  src_vec -> ExtractView(&src_vec_ptr);
  MPI_Allreduce(MPI_IN_PLACE, src_vec_ptr, num_ditches_, MPI_DOUBLE, MPI_SUM, comm);


  //std::cout <<"After MPI_Allreduce\n";
  //std::cout << *src_vec<<"\n";

  //double norm;

  //std::cout<<sub_sink<<"\n";
  // sub_sink.Norm2(&norm);
  // if (norm >1e-10) {
  //   //sub_sink.Norm2(&norm);
  //   std::cout << "Sink norm "<<norm <<"\n";
  //   //    exit(0);
  // }

}

// Required methods from SecondaryVariableFieldEvaluator
void
DistributedTilesRateEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
                     Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{

  result->PutScalar(0.0);
  
}

void
DistributedTilesRateEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S) {

  //sources_key_ = Key("subdomain_sources");
  if (!S->HasField(sources_key_)){
    S->RequireConstantVector(sources_key_, num_ditches_ * num_component_);
    Teuchos::RCP<Field> field =  S->GetField(sources_key_, "state");
    Teuchos::RCP<Field_ConstantVector> cvfield =
      Teuchos::rcp_dynamic_cast<Field_ConstantVector>(field, true);
    cvfield->CreateData();

    field->GetConstantVectorData()->PutScalar(0.0);
    field->set_initialized();

    if (!implicit_) {
      Key prev_step("prev_timestep");
      field->RequireCopy(prev_step);
    }
  }
  
  SecondaryVariableFieldEvaluator::EnsureCompatibility(S);
  S->GetField(my_key_, my_key_)->set_io_vis();
}

void
DistributedTilesRateEvaluator::InitializeFromPlist_() {

  
  
}



} //namespace
} //namespace
} //namespace

