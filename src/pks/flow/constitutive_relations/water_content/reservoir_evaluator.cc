/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/
//! Incorporates fluxes out of a reservoir/water body.
/*!


.. _field-evaluator-type-reservoir-spec:
.. admonition:: field-evaluator-type-reservoir-spec

   DEPENDENCIES:

   - `"water content`"
   - `"cell volume`"

*/

#include "reservoir_evaluator.hh"
#include "reservoir_model.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

// Constructor from ParameterList
ReservoirEvaluator::ReservoirEvaluator(Teuchos::ParameterList& plist) :
    EvaluatorSecondaryMonotypeCV(plist)
{
  InitializeFromPlist_();
}


// Virtual copy constructor
Teuchos::RCP<Evaluator>
ReservoirEvaluator::Clone() const
{
  return Teuchos::rcp(new ReservoirEvaluator(*this));
}


// Initialize by setting up dependencies
void
ReservoirEvaluator::InitializeFromPlist_()
{
  // require regions
  region_waterbody_ = plist_.get<std::string>("reservoir water body region");
  region_outlet_ = plist_.get<std::string>("reservoir outlet region");

  // Set up my dependencies
  // - defaults to prefixed via domain
  Key domain_name = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  // - pull Keys from plist
  // dependency: porosity
  wc_key_ = Keys::readKey(plist_, domain_name, "water content", "water_content");
  dependencies_.insert(KeyTag{wc_key_, tag});

  // dependency: cell_volume
  cv_key_ = Keys::readKey(plist_, domain_name, "cell volume", "cell_volume");
  dependencies_.insert(KeyTag{cv_key_, tag});
}


void
ReservoirEvaluator::Evaluate_(const State& S,
        const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  const auto& wc = *S.Get<CompositeVector>(wc_key_, tag).ViewComponent("cell", false);
  const auto& cv = *S.Get<CompositeVector>(cv_key_, tag).ViewComponent("cell", false);
  auto& result_v = *result[0]->ViewComponent("cell",false);

  AmanziMesh::Entity_ID_List waterbody_cells, outlet_cells;
  result[0]->Mesh()->get_set_entities(region_waterbody_, AmanziMesh::Entity_kind::CELL,
          AmanziMesh::Parallel_type::OWNED, &waterbody_cells);
  result[0]->Mesh()->get_set_entities(region_outlet_, AmanziMesh::Entity_kind::CELL,
          AmanziMesh::Parallel_type::OWNED, &outlet_cells);

  // compute the local water body quantities
  double input_local[2] = { 0., 0. };
  for (AmanziMesh::Entity_ID c : waterbody_cells) {
    input_local[0] += wc[0][c];
    input_local[1] += cv[0][c];
  }

  double input_global[2];
  result[0]->Mesh()->get_comm()->SumAll(input_local, input_global, 2);

  // if there is no water to give, cut out now
  if (input_global[0] < 1.e-5) return;

  // compute on the outlet node
  double output_local[2] = { 0., 0. };
  if (outlet_cells.size() > 0) {
    if (model_ == Teuchos::null) {
      Teuchos::ParameterList& sublist = plist_.sublist("reservoir model parameters");
      model_ = createReservoirModel(sublist);
    }

    // expected in units of mols/s
    output_local[1] = outlet_cells.size();
    output_local[0] = model_->computeDischarge(input_global[0]);

    // set the source, units of mols/m^2/s
    result_v[0][outlet_cells[0]] = output_local[0] / cv[0][outlet_cells[0]];
  }

  // set the sink -- must communicate
  double output_global[2] = { 0., 0. };
  result[0]->Mesh()->get_comm()->SumAll(output_local, output_global, 2);
  if (output_global[1] != 1) {
    Errors::Message msg;
    msg << "ReservoirEvaluator: model expects that only one cell is the outlet of the reservoir, but outlet region \"" << region_outlet_ << "\" does not satisfy that criteria.";
    Exceptions::amanzi_throw(msg);
  }

  // weight the source by water content, convert units to mol/m^2/s, and
  // negative due to sink
  for (AmanziMesh::Entity_ID c : waterbody_cells) {
    result_v[0][c] = -output_global[0] * wc[0][c] / input_global[0] / cv[0][c];
  }
}


} //namespace
} //namespace
} //namespace
