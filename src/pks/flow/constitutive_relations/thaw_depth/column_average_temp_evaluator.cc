/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The column average temperature evaluator gets the subsurface temperature and number of cells (related to depth), and returns the average column temperature.

  Authors: Ahmad Jan (jana@ornl.gov)
*/

#include "column_average_temp_evaluator.hh"

namespace Amanzi {
namespace Flow {



ColumnAverageTempEvaluator::ColumnAverageTempEvaluator(Teuchos::ParameterList& plist)
    : SecondaryVariableFieldEvaluator(plist)
{
  domain_ = Keys::getDomain(my_key_);
  auto pos = domain_.find_last_of('_');
  int col_id = std::stoi(domain_.substr(pos+1, domain_.size()));
  
  std::stringstream domain_ss;
  domain_ss << "column_"<< col_id;
  temp_key_ = Keys::getKey(domain_ss.str(),"temperature");
  dependencies_.insert(temp_key_);
  
  depth_ = plist_.get<double>("depth from surface [m]", 0); // depth from the surface
  ncells_depth_ = plist_.get<int>("number of cells [m]", -1); // or number of cells
}
  

Teuchos::RCP<FieldEvaluator>
ColumnAverageTempEvaluator::Clone() const
{
  return Teuchos::rcp(new ColumnAverageTempEvaluator(*this));
}


void
ColumnAverageTempEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{ 
  Epetra_MultiVector& res_c = *result->ViewComponent("cell",false);
  
  // search through the column and find the deepest unfrozen cell

  std::string domain_ss = Keys::getDomain(temp_key_);
  const auto& top_z_centroid = S->GetMesh(domain_ss)->face_centroid(0);
  AmanziGeometry::Point z_centroid(top_z_centroid);

  const auto& temp_c = *S->GetFieldData(temp_key_)->ViewComponent("cell", false);
  
  int col_cells = temp_c.MyLength();
  double temp_sum = 0;
  int count = 0 ;
  
  AMANZI_ASSERT (ncells_depth_ <= col_cells);

  for (int i=0; i!=col_cells; ++i) {
    if (depth_ > 0.0) {
      z_centroid = S->GetMesh(domain_ss)->face_centroid(i+1);
      double z_depth = top_z_centroid[2] - z_centroid[2];
      if (z_depth <= depth_) {
        temp_sum += temp_c[0][i];
        count += 1;
      }
      else {
        break;
      }
    }
    else if (ncells_depth_ > 0 && i <= ncells_depth_) {
      temp_sum += temp_c[0][i];
      count += 1;
    }
  }
  
  res_c[0][0] = count >0 ? temp_sum/count : 0.0;

}
  
void
ColumnAverageTempEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
               Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{}

 
// Custom EnsureCompatibility forces this to be updated once.
bool
ColumnAverageTempEvaluator::HasFieldChanged(const Teuchos::Ptr<State>& S,
        Key request)
{
  bool changed = SecondaryVariableFieldEvaluator::HasFieldChanged(S,request);

  if (!updated_once_) {
    UpdateField_(S);
    updated_once_ = true;
    return true;
  }
  return changed;
}

void
ColumnAverageTempEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S)
{

  AMANZI_ASSERT(my_key_ != std::string(""));
   
  Teuchos::RCP<CompositeVectorSpace> my_fac = S->RequireField(my_key_, my_key_);
  
  // check plist for vis or checkpointing control
  bool io_my_key = plist_.get<bool>(std::string("visualize ")+my_key_, true);
  S->GetField(my_key_, my_key_)->set_io_vis(io_my_key);
  bool checkpoint_my_key = plist_.get<bool>(std::string("checkpoint ")+my_key_, false);
  S->GetField(my_key_, my_key_)->set_io_checkpoint(checkpoint_my_key);
  
  if (my_fac->Mesh() != Teuchos::null) {
    // Recurse into the tree to propagate info to leaves.
    for (KeySet::const_iterator key=dependencies_.begin();
         key!=dependencies_.end(); ++key) {
      S->RequireFieldEvaluator(*key)->EnsureCompatibility(S);
    }
  }
}


} //namespace
} //namespace
