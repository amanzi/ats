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
#include "overland_mass_source_evaluator.hh"
#include "FunctionFactory.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {


OverlandMassSourceEvaluator::OverlandMassSourceEvaluator(Teuchos::ParameterList& plist) : EvaluatorSecondaryMonotypeCV(plist)
{
  domain_ = Keys::getDomain(my_keys_.front().first);
  auto tag = my_keys_.front().second;

  cv_key_ = Keys::readKey(plist, domain_, "cell volume", "cell_volume");
  dependencies_.insert(KeyTag{ cv_key_, tag });
  molar_density_key_ = Keys::readKey(plist, domain_, "molar density liquid", "molar_density_liquid");
  dependencies_.insert(KeyTag{ molar_density_key_, tag });
  field_src_key_ = Keys::readKey(plist, domain_, "overland mass source", "overland_mass_source");
  // dependencies_.insert(KeyTag{ field_src_key_, tag });

  Teuchos::ParameterList& source_func = plist.sublist("function");
  FunctionFactory fac;
  QC_curve_ = Teuchos::rcp(fac.Create(source_func));

}

// Required methods from SecondaryVariableFieldEvaluator

void
OverlandMassSourceEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result) 
{

  Tag tag = my_keys_.front().second;
  
  const auto& cv = *S.Get<CompositeVector>(cv_key_, tag).ViewComponent("cell", false);
  const auto& molar_den = *S.Get<CompositeVector>(molar_density_key_, tag).ViewComponent("cell", false);
  const auto& water_from_field = *S.Get<CompositeVector>(field_src_key_, tag).ViewComponent("face", false);  
  auto& surf_src= *result[0]->ViewComponent("cell"); // not being reference

  const AmanziMesh::Mesh& mesh = *result[0]->Mesh();
  std::unordered_map<int, std::vector<int>> face_to_cells;
  AmanziMesh::Entity_ID ncells = cv.MyLength();

  // Populate the face_to_cells map
  for (AmanziMesh::Entity_ID c = 0; c != ncells; ++c) {
    const auto& [faces, dirs] = mesh.getCellFacesAndDirections(c);
    for (int f : faces) {
        face_to_cells[f].push_back(c);
    }    
  }

  // Output for each cell which faces are shared or not shared
  for (AmanziMesh::Entity_ID c = 0; c != ncells; ++c) {
    double total_external_flux = 0;
    const auto& [faces, dirs] = mesh.getCellFacesAndDirections(c);    
    for (int f : faces) {
      if (face_to_cells[f].size() == 1) { //face is not shared, hence external
        total_external_flux += water_from_field[0][f];
      }
    }
    // convert from mol/s to m3/s
    double total_flux_meter = total_external_flux * cv[0][c] / molar_den[0][c];

    // transport source (mass) as a function of discharge (e.g. overland)
    double source_mass = (*QC_curve_)(std::vector<double>{total_flux_meter});
    surf_src[0][c] = source_mass;
  }
}

} // namespace Relations
} // namespace Flow
} // namespace Amanzi 