/* -*-  mode: c++; indent-tabs-mode: nil -*- */
//! SubgridAggregateEvaluator restricts a field to the subgrid version of the same field.

/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/


/*!

 * `"source domain name`" ``[string]`` Domain name of the source mesh.

ONE OF:
* `"field suffix`" ``[string]`` Set the suffix of the variable
OR
* `"field key`" ``[string]`` **DOMAIN-FIELD_SUFFIX** 

  
 */

#include "SubgridAggregateEvaluator.hh"

namespace Amanzi {
namespace Relations {

SubgridAggregateEvaluator::SubgridAggregateEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist)
{
  domain_ = Keys::getDomain(my_key_);
  source_domain_ = plist_.get<std::string>("source domain name");
  var_key_ = Keys::getVarName(my_key_);
  
}

Teuchos::RCP<FieldEvaluator>
SubgridAggregateEvaluator::Clone() const
{
  return Teuchos::rcp(new SubgridAggregateEvaluator(*this));
}

// Required methods from SecondaryVariableFieldEvaluator
void
SubgridAggregateEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{

  int ncells = S->GetMesh("surface_star")->num_entities(AmanziMesh::CELL,
                                                        AmanziMesh::Parallel_type::OWNED);

  for (int c =0; c < ncells; c++) {
    std::stringstream name;
    int id = S->GetMesh(domain_)->cell_map(false).GID(c);
    name << "surface_column_"<< id;
    Key source_key = Keys::getKey(name.str(),var_key_);

    const auto& source = *S->GetFieldData(source_key)->ViewComponent("cell",false);
    
    AMANZI_ASSERT(source.MyLength() == 1);
    (*result->ViewComponent("cell", false))[0][c] = source[0][0];
    
  }

}

void
SubgridAggregateEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  result->PutScalar(1.);
}


void
SubgridAggregateEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S)
{ 

  int ncells = S->GetMesh("surface_star")->num_entities(AmanziMesh::CELL,
                                                        AmanziMesh::Parallel_type::OWNED);
  
  for (int c =0; c < ncells; c++){
    std::stringstream name;
    int id = S->GetMesh("surface_star")->cell_map(false).GID(c);
    name << "surface_column_"<< id;
    Key temp_key = Keys::getKey(name.str(),var_key_);
    dependencies_.insert(temp_key);
  }

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
      S->RequireField(*key)->SetMesh(S->GetMesh(Keys::getDomain(*key)))
        ->AddComponent("cell", AmanziMesh::CELL, 1);
    S->RequireFieldEvaluator(*key)->EnsureCompatibility(S);
    }
  }

}


} // namespace
} // namespace

