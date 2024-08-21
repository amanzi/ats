/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Phong V.V. Le (lepv@ornl.gov)
*/

#include "Key.hh"
#include "Factory.hh"
#include "Function.hh"
#include "qc_relation_overland_evaluator.hh"
#include "FunctionFactory.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {


QCRelationOverlandEvaluator::QCRelationOverlandEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  domain_ = Keys::getDomain(my_keys_.front().first);
  auto tag = my_keys_.front().second;

  cv_key_ = Keys::readKey(plist, domain_, "cell volume", "cell_volume");
  dependencies_.insert(KeyTag{ cv_key_, tag });
  molar_density_key_ =
    Keys::readKey(plist, domain_, "molar density liquid", "molar_density_liquid");
  dependencies_.insert(KeyTag{ molar_density_key_, tag });
  field_src_key_ = Keys::readKey(plist, domain_, "overland source", "water_flux");

  // Create a Q-C curve using "function" in parameter list
  Teuchos::ParameterList& source_func = plist.sublist("function");
  FunctionFactory factory;
  QC_curve_ = Teuchos::rcp(factory.Create(source_func));
}

// Required methods from SecondaryVariableFieldEvaluator
void
QCRelationOverlandEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;

  const auto& cv = *S.Get<CompositeVector>(cv_key_, tag).ViewComponent("cell", false);
  const auto& molar_den =
    *S.Get<CompositeVector>(molar_density_key_, tag).ViewComponent("cell", false);
  const auto& water_from_field =
    *S.Get<CompositeVector>(field_src_key_, tag).ViewComponent("face", false);
  auto& surf_src = *result[0]->ViewComponent("cell"); // not being reference

  const AmanziMesh::Mesh& mesh = *result[0]->Mesh();
  AmanziMesh::Entity_ID ncells = cv.MyLength();

  // We consider flux only from hillslope to river channel cells. This flux brings solute source to river
  // We separate: (i) external faces which connect river cell and hill slope cell and (ii) internal faces which connect two river cells
  // mesh.getFaceCells(f).size() == 1: Face on the river bank connecting overland and river cells (use)
  // mesh.getFaceCells(f).size() == 2: Internal face connecting two river cells (ignore)
  for (AmanziMesh::Entity_ID c = 0; c != ncells; ++c) {
    double total_external_flux = 0;
    const auto& [faces, dirs] = mesh.getCellFacesAndDirections(c);
    int nfaces = faces.size();
    
    for (int i = 0; i < nfaces; i++) {
      int f = faces[i];     // get each face of the cell
      double dir = dirs[i]; // 1: water goes out of the cell; -1: water goes into the cell

      if (mesh.getFaceCells(f).size() == 1) {
        // external faces which contributes sinks/sources
        total_external_flux += std::abs(water_from_field[0][f]) * (-dir);
      }
    }
    
    // convert from mol/s to m3/s. We do NOT multiply with cv here
    double total_flux_meter = total_external_flux / molar_den[0][c];

    // transport source (mass) as a function of discharge (e.g. overland)
    double source_transport = (*QC_curve_)(std::vector<double>{ total_flux_meter });
    surf_src[0][c] = source_transport * total_flux_meter;
  }
}

} // namespace Relations
} // namespace Flow
} // namespace Amanzi