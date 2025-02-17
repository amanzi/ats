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
  const auto& water_from_field =
    *S.Get<CompositeVector>(field_src_key_, tag).ViewComponent("cell", false);
  auto& surf_src = *result[0]->ViewComponent("cell"); // not being reference
  const AmanziMesh::Mesh& mesh = *result[0]->Mesh();
  double field_flow, source_mass;
  
  // Loop through each cell
  AmanziMesh::Entity_ID ncells = cv.MyLength();
  for (AmanziMesh::Entity_ID c = 0; c != ncells; ++c) {
    const auto& [faces, dirs] = mesh.getCellFacesAndDirections(c);
    int nfaces = faces.size();
    if (nfaces < 4) {   // skip triangle cells that are not river corridor
      surf_src[0][c] = 0;
    } else {
      if (extensive_) {
        // convert extensive quantity in mol/m2/s to m3/s
        field_flow = water_from_field[0][c] / molar_den[0][c];
      } else {
        // convert intensive quantity from mol/s to m3/s
        field_flow = water_from_field[0][c] * cv[0][c] / molar_den[0][c];
      }

      // transport source (concentration g/m3) as a function of discharge from a field (e.g. tile, groundwater)
      source_mass = (*QC_curve_)(std::vector<double>{field_flow});

      // return solute mass rate by multiplying with discharge (molC/s)
      // Here we assume the molar mass is 1. TODO: add molar mass to the function.
      surf_src[0][c] = source_mass * field_flow / cv[0][c];
    }
  }
}

} // namespace Relations
} // namespace Flow
} // namespace Amanzi 