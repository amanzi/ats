/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Phong V.V. Le (lepv@ornl.gov)
*/
#include <algorithm>
#include "Key.hh"
#include "Factory.hh"
#include "Function.hh"
#include "qc_relation_field_evaluator.hh"
#include "FunctionFactory.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {


QCRelationFieldEvaluator::QCRelationFieldEvaluator(Teuchos::ParameterList& plist) : EvaluatorSecondaryMonotypeCV(plist)
{
  domain_ = Keys::getDomain(my_keys_.front().first);
  auto tag = my_keys_.front().second;

  cv_key_ = Keys::readKey(plist, domain_, "cell volume", "cell_volume");
  dependencies_.insert(KeyTag{ cv_key_, tag });
  molar_density_key_ = Keys::readKey(plist, domain_, "molar density liquid", "molar_density_liquid");
  tcc_key_ = Keys::readKey(plist, domain_, "concentration", "total_component_concentration");
  dependencies_.insert(KeyTag{ molar_density_key_, tag });
  field_src_key_ = Keys::readKey(plist, domain_, "field source", "water_source_field");
  dependencies_.insert(KeyTag{ field_src_key_, tag });

  // Check extensive or intensive quantity
  extensive_ = plist.get<bool>("extensive", false);

  // Create a Q-C curve using "function" in parameter list
  Teuchos::ParameterList& qc_list_func = plist.sublist("function");
  FunctionFactory factory;
  QC_curve_ = Teuchos::rcp(factory.Create(qc_list_func));
}

// Required methods from SecondaryVariableFieldEvaluator
void
QCRelationFieldEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result) 
{

  Tag tag = my_keys_.front().second;

  const auto& cv = *S.Get<CompositeVector>(cv_key_, tag).ViewComponent("cell", false);
  const auto& molar_den =
    *S.Get<CompositeVector>(molar_density_key_, tag).ViewComponent("cell", false);
  const auto& tcc = *S.Get<CompositeVector>(tcc_key_, tag).ViewComponent("cell", false);    
  const auto& water_from_field =
    *S.Get<CompositeVector>(field_src_key_, tag).ViewComponent("cell", false);
  auto& surf_src = *result[0]->ViewComponent("cell"); // not being reference
  const AmanziMesh::Mesh& mesh = *result[0]->Mesh();
  double field_flow, source_mass, tcc_current;
  
  // Loop through each cell
  AmanziMesh::Entity_ID ncells = cv.MyLength();
  for (AmanziMesh::Entity_ID c = 0; c != ncells; ++c) {
    
    if (extensive_) {
      // convert extensive quantity in mol/s to m3/s
      field_flow = water_from_field[0][c] / molar_den[0][c];
    } else {
      // convert intensive quantity from mol/m2/s to m3/s
      field_flow = water_from_field[0][c] * cv[0][c] / molar_den[0][c];
    }
    // current concentration
    tcc_current = tcc[0][c];

    // transport source (concentration g/m3) as a function of discharge from a field (e.g. tile, groundwater)
    source_mass = (*QC_curve_)(std::vector<double>{std::abs(field_flow)});
    
    // temporarily assume molar mass is 1.0
    // TODO: read molar mass from the xml file
    if (field_flow > 0) {
      // positive flux means source, concentration is from source_transport
      surf_src[0][c] = source_mass * field_flow * 1.0;
    } else {
      // negative flux means sink, concentration is the same as the current concentration
      surf_src[0][c] = tcc_current * (-field_flow) * 1.0;
    }    
  }
}

} // namespace Relations
} // namespace Flow
} // namespace Amanzi 