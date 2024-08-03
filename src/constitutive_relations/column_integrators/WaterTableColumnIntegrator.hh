/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*

This is an evaluator for integrating arbitrary functionals across columns to
get a related quantity.  Example uses might be computing the column-averaged
temperature, finding the depth to water table, or similar.

Clients should provide a struct functor that does the actual work, and returns
-1 if the loop over columns should break.

*/

#pragma once

#include "Factory.hh"
#include "EvaluatorColumnIntegrator.hh"
#include <math.h>

namespace Amanzi {
namespace Relations {

template <class Parser, class Integrator>
class WaterTableColumnIntegrator : public EvaluatorColumnIntegrator<Parser, Integrator> {
 public:
  using EvaluatorColumnIntegrator<Parser, Integrator>::EvaluatorColumnIntegrator;
  WaterTableColumnIntegrator(const WaterTableColumnIntegrator& other) = default;
  Teuchos::RCP<Evaluator> Clone() const override;

  // Required methods from EvaluatorColumnIntegrator to overide
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;

  Teuchos::ParameterList plist_ = EvaluatorColumnIntegrator<Parser, Integrator>::plist_;
  KeyTagSet dependencies_ = EvaluatorColumnIntegrator<Parser, Integrator>::dependencies_;

 private:
  static Utils::RegisteredFactory<Evaluator, WaterTableColumnIntegrator<Parser, Integrator>> reg_;
};


template <class Parser, class Integrator>
Teuchos::RCP<Evaluator>
WaterTableColumnIntegrator<Parser, Integrator>::Clone() const
{
  return Teuchos::rcp(new WaterTableColumnIntegrator<Parser, Integrator>(*this));
}


// Required methods from EvaluatorColumnIntegrator
template <class Parser, class Integrator>
void
WaterTableColumnIntegrator<Parser, Integrator>::Evaluate_(
  const State& S,
  const std::vector<CompositeVector*>& result)
{
  // collect the dependencies and mesh, and instantiate the integrator functor
  std::vector<const Epetra_MultiVector*> deps;
  for (const auto& dep : dependencies_) {
    deps.emplace_back(
      S.Get<CompositeVector>(dep.first, dep.second).ViewComponent("cell", false).get());
  }
  auto mesh = result[0]->Mesh()->getParentMesh();
  Integrator integrator(plist_, deps, &*mesh);

  Epetra_MultiVector& res = *result[0]->ViewComponent("cell", false);

  for (int col = 0; col != res.MyLength(); ++col) {
    // for each column, loop over cells calling the integrator until stop is
    // requested or the column is complete
    AmanziGeometry::Point val(0., 0., NAN);
    auto col_cell = mesh->columns.getCells(col);
    double h_top = mesh->getCellCentroid(col_cell[0])[2]
        + mesh->getCellVolume(col_cell[0]) * integrator.coefficient(col) / 2;
    double h_bot = mesh->getCellCentroid(col_cell[col_cell.size() - 1])[2]
        - mesh->getCellVolume(col_cell[col_cell.size() - 1]) * integrator.coefficient(col) / 2;
    double h_end, h_half0, h_half1;
    h_end = mesh->getCellCentroid(col_cell[0])[2]; // default at top centroid
    h_half0 = mesh->getCellVolume(col_cell[col_cell.size() - 1]) * integrator.coefficient(col) / 2;
    h_half1 = mesh->getCellVolume(col_cell[0]) * integrator.coefficient(col) / 2 * (-1);
    for (int i = col_cell.size() - 1; i >= 0; --i) { // loop from bottom up looking for the 1st unsaturated cell 
      bool completed = integrator.scan(col, col_cell[i], val);
      if (completed) {
        h_end = mesh->getCellCentroid(col_cell[i])[2]; // the first unsaturated cell centroid from bottom up
        break;
      }
    }

    // Use val[2] to track centroid, and val[1], val[0] to track cell pressure 
    // or volume determined by using interpolation or not. 
    if (std::isnan(val[2])) { // completed at first loop cell
      res[0][col] = h_top - h_end + h_half0;
    } else if (val[2] == h_end) { // fail to find satisfied cell util end of loop 
      res[0][col] = h_top - val[2] + h_half1;
    } else {
      if (plist_.get<bool>("determined by pressure interpolation")) {
        res[0][col] = (val[2] - h_end) * (101325. - val[0]) / (val[1] - val[0]) 
                    + (h_top - val[2]);
      } else {
        res[0][col] = h_top - h_bot - integrator.coefficient(col) * val[0];
      }
    }
  }
}

} // namespace Relations
} // namespace Amanzi
