/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Phong V.V. Le
*/

#include "Key.hh"
#include "Factory.hh"
#include "Function.hh"
#include "field_mass_source_evaluator.hh"
#include "FunctionFactory.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {


FieldMassSourceEvaluator::FieldMassSourceEvaluator(Teuchos::ParameterList& plist) : EvaluatorSecondaryMonotypeCV(plist)
{
  domain_ = Keys::getDomain(my_keys_.front().first);
  auto tag = my_keys_.front().second;

  cv_key_ = Keys::readKey(plist, domain_, "cell volume", "cell_volume");
  dependencies_.insert(KeyTag{ cv_key_, tag });
  molar_density_key_ = Keys::readKey(plist, domain_, "molar density liquid", "molar_density_liquid");
  dependencies_.insert(KeyTag{ molar_density_key_, tag });
  field_src_key_ = Keys::readKey(plist, domain_, "field mass source", "field_mass_source");
  dependencies_.insert(KeyTag{ field_src_key_, tag });

  Teuchos::ParameterList& gate_func = plist.sublist("function");
  FunctionFactory fac;
  QC_curve_ = Teuchos::rcp(fac.Create(gate_func));

}

// Required methods from SecondaryVariableFieldEvaluator

void
FieldMassSourceEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result) 
{

  Tag tag = my_keys_.front().second;
  
  const auto& cv = *S.Get<CompositeVector>(cv_key_, tag).ViewComponent("cell", false);
  const auto& molar_den = *S.Get<CompositeVector>(molar_density_key_, tag).ViewComponent("cell", false);
  const auto& water_from_field = *S.Get<CompositeVector>(field_src_key_, tag).ViewComponent("cell", false);
  
  auto& surf_src= *result[0]->ViewComponent("cell"); // not being reference

  double total = 0.0;
  const AmanziMesh::Mesh& mesh = *result[0]->Mesh();

  AmanziMesh::Entity_ID ncells = cv.MyLength();
  for (AmanziMesh::Entity_ID c = 0; c != ncells; ++c) {
    // convert discharge from mol/s to m3/s
    double field_flow = water_from_field[0][c] * cv[0][c] / molar_den[0][c];

    // transport source (mass) as a function of particular discharge (e.g. tile, overland)
    double source_mass = (*QC_curve_)(std::vector<double>{field_flow});
    surf_src[0][c] = source_mass;
  }
}

} // namespace Relations
} // namespace Flow
} // namespace Amanzi 