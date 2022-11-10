/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*

  Evaluates water/solute sinks which represent distributed tiles in subsurface 
  
  License: see $ATS_DIR/COPYRIGHT
  Author: Daniil Svyatsky (dasvyat@lanl.gov) 
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
    EvaluatorSecondaryMonotypeCV(plist),
    compatibility_checked_(false)
{

  domain_ = Keys::getDomain(my_keys_.front().first);
  auto tag = my_keys_.front().second;

  subsurface_marks_key_ = Keys::readKey(plist, domain_, "catchments_id", "catchments_id");
  dist_sources_key_ = plist.get<std::string>("accumulated source key", "subdomain_sources");
  factor_key_ = plist.get<std::string>("factor field key", "");
  num_component_ = plist.get<int>("number of components", 1);
    
  
  num_ditches_ = plist.get<int>("number of ditches");
  p_enter_ = plist.get<double>("entering pressure", 101325);
  k_ = plist.get<double>("tile permeability");
  implicit_ = plist.get<bool>("implicit drainage", true);

  dependencies_.insert(KeyTag{subsurface_marks_key_, tag});
  pres_key_ = Keys::readKey(plist, domain_, "pressure", "pressure");
  dependencies_.insert(KeyTag{pres_key_, tag});
  mol_dens_key_ = Keys::readKey(plist, domain_, "molar density", "molar_density_liquid");
  dependencies_.insert(KeyTag{mol_dens_key_, tag});

  if (factor_key_ != "") dependencies_.insert(KeyTag{factor_key_, tag});
  
}

// Required methods from SecondaryVariableFieldEvaluator
void
DistributedTilesRateEvaluator::Evaluate_(const State& S,
                                         const std::vector<CompositeVector*>& result)
{

  auto tag = my_keys_.front().second;
  // double t0 = S->time();
  double dt = S.Get<double>("dt", tag);
  // //  double t1 = S->final_time();


  const auto& pres = *S.GetPtr<CompositeVector>(pres_key_, tag)->ViewComponent("cell", false);
  const auto& dens = *S.GetPtr<CompositeVector>(mol_dens_key_, tag)->ViewComponent("cell", false);
  
  const auto& cv =
    *S.GetPtr<CompositeVector>(Keys::getKey(domain_,"cell_volume"), tag)->ViewComponent("cell",false);
  const auto& sub_marks = *S.GetPtr<CompositeVector>(subsurface_marks_key_, tag)->ViewComponent("cell", false);
  
  auto& sub_sink = *result[0]->ViewComponent("cell",false);
  sub_sink.PutScalar(0);
  
  AmanziMesh::Entity_ID ncells = sub_marks.MyLength();

  double total = 0.0;
  dist_src_vec_->PutScalar(0.0);
  sub_sink.PutScalar(0.0);
  int num_vectors = 1;
  int test = sub_sink.NumVectors();

  if (factor_key_!=""){
    num_vectors = S.GetPtr<CompositeVector>(factor_key_, tag)->ViewComponent("cell",false)->NumVectors();
    AMANZI_ASSERT(num_vectors == sub_sink.NumVectors());
    AMANZI_ASSERT(num_vectors == num_component_);
  }

  if (abs(dt) > 1e-13) {
    const Epetra_MultiVector *factor = nullptr;
    if (factor_key_!=""){
      factor = &(*(S.GetPtr<CompositeVector>(factor_key_, tag)->ViewComponent("cell", false)));
    }
    for (AmanziMesh::Entity_ID c=0; c!=ncells; ++c) {
      if (sub_marks[0][c] > 0) {
        double val = std::min(p_enter_ - pres[0][c], 0.0)*dens[0][c]*k_;
        //val = 1e-6;
        //std::cout<<"sink val "<<" "<<val<<" "<< pres[0][c]<<"\n";

        if (factor_key_!=""){
          for (int i=0; i<num_component_; ++i){
            sub_sink[i][c] = (*factor)[i][c] * val;
            (*dist_src_vec_)[sub_marks[0][c] - 1 + i* num_ditches_] += (*factor)[i][c] * val * dt * cv[0][c];
            total = total + (*factor)[i][c] * val * dt * cv[0][c];
          }
        }else{
          sub_sink[0][c] = val;
          (*dist_src_vec_)[sub_marks[0][c] - 1] +=  val * dt * cv[0][c];
          total = total + val * dt * cv[0][c];
        }
      }
    }
  }
  // std::cout <<"Sink vector\n";
  // std::cout << *dist_src_vec_<<"\n";
  //std::cout<<"Total sink "<<my_key_<<" "<<total<<"\n";
  // total = 0.;
  // for (AmanziMesh::Entity_ID c=0; c!=ncells; ++c) {
  //   total += sub_sink[0][c] * cv[0][c];
  // }
  // std::cout<<"Total sink field "<<" "<<total<<"\n";

  
  Teuchos::RCP<const Comm_type> comm_p = S.GetMesh(domain_)->get_comm();
  Teuchos::RCP<const MpiComm_type> mpi_comm_p =
    Teuchos::rcp_dynamic_cast<const MpiComm_type>(comm_p);
  const MPI_Comm& comm = mpi_comm_p->Comm();

  double *src_vec_ptr;
  dist_src_vec_ -> ExtractView(&src_vec_ptr);
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
DistributedTilesRateEvaluator::EvaluatePartialDerivative_(const State& S,
    const Key& wrt_key, const Tag& wrt_tag, const std::vector<CompositeVector*>& result)
{
  result[0]->PutScalar(0.0);
}

void
DistributedTilesRateEvaluator::EnsureCompatibility_Structure_(State& S) {

  auto tag = my_keys_.front().second;
  if (!S.HasRecord(dist_sources_key_, tag)){
    S.Require<Epetra_Vector, Epetra_Vector_Factory>(dist_sources_key_, tag, "state").set_size(num_ditches_);
    S.GetRecordSetW(dist_sources_key_).CreateData();    
  }

  dist_src_vec_ = S.GetPtrW<Epetra_Vector>(dist_sources_key_, tag, "state");
  dist_src_vec_->PutScalar(0.0);
  S.GetRecordW(dist_sources_key_, tag, "state").set_initialized();
  
  EvaluatorSecondaryMonotypeCV::EnsureCompatibility_Structure_(S);
  
}

void
DistributedTilesRateEvaluator::InitializeFromPlist_() {

  
  
}



} //namespace
} //namespace
} //namespace

