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
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Relations {

template <class Parser, class Integrator>
class EvaluatorColumnIntegrator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit EvaluatorColumnIntegrator(Teuchos::ParameterList& plist);
  EvaluatorColumnIntegrator(const EvaluatorColumnIntegrator& other) = default;
  Teuchos::RCP<Evaluator> Clone() const override;

  // Disables derivatives
  virtual bool
  IsDifferentiableWRT(const State& S, const Key& wrt_key, const Tag& wrt_tag) const override;

 protected:
  // Implements custom EC to use dependencies from subsurface for surface
  // vector.
  virtual void EnsureCompatibility_ToDeps_(State& S) override;

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

 private:
  static Utils::RegisteredFactory<Evaluator, EvaluatorColumnIntegrator<Parser, Integrator>> reg_;
};


template <class Parser, class Integrator>
EvaluatorColumnIntegrator<Parser, Integrator>::EvaluatorColumnIntegrator(
  Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  Parser parser(plist_, my_keys_.front());
  dependencies_ = std::move(parser.dependencies);
}


template <class Parser, class Integrator>
Teuchos::RCP<Evaluator>
EvaluatorColumnIntegrator<Parser, Integrator>::Clone() const
{
  return Teuchos::rcp(new EvaluatorColumnIntegrator<Parser, Integrator>(*this));
}


// Implements custom EC to use dependencies from subsurface for surface
// vector.
template <class Parser, class Integrator>
void
EvaluatorColumnIntegrator<Parser, Integrator>::EnsureCompatibility_ToDeps_(State& S)
{
  const auto& fac = S.Require<CompositeVector, CompositeVectorSpace>(my_keys_.front().first,
                                                                     my_keys_.front().second);
  if (fac.Mesh() != Teuchos::null) {
    CompositeVectorSpace dep_fac;
    dep_fac.SetMesh(fac.Mesh()->getParentMesh())
      ->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

    for (const auto& dep : dependencies_) {
      if (Keys::getDomain(dep.first) == Keys::getDomain(my_keys_.front().first)) {
        S.Require<CompositeVector, CompositeVectorSpace>(dep.first, dep.second).Update(fac);
      } else {
        S.Require<CompositeVector, CompositeVectorSpace>(dep.first, dep.second).Update(dep_fac);
      }
    }
  }
}


// Disables derivatives
template <class Parser, class Integrator>
bool
EvaluatorColumnIntegrator<Parser, Integrator>::IsDifferentiableWRT(const State& S,
                                                                   const Key& wrt_key,
                                                                   const Tag& wrt_tag) const
{
  return false;
}


// Required methods from EvaluatorSecondaryMonotypeCV
template <class Parser, class Integrator>
void
EvaluatorColumnIntegrator<Parser, Integrator>::Evaluate_(
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
    AmanziGeometry::Point val(0., 0.);
    auto col_cell = mesh->columns.getCells(col);
    for (int i = 0; i != col_cell.size(); ++i) {
      bool completed = integrator.scan(col, col_cell[i], val);
      if (completed) break;
    }

    // val[1] is typically e.g. cell volume, but can be 0 to indicate no
    // denominator.  Coefficient provides a hook for column-wide multiples
    // (e.g. 1/surface area).
    if (val[1] > 0.)
      res[0][col] = integrator.coefficient(col) * val[0] / val[1];
    else
      res[0][col] = integrator.coefficient(col) * val[0];
  }
}


// Required methods from EvaluatorSecondaryMonotypeCV
template <class Parser, class Integrator>
void
EvaluatorColumnIntegrator<Parser, Integrator>::EvaluatePartialDerivative_(
  const State& S,
  const Key& wrt_key,
  const Tag& wrt_tag,
  const std::vector<CompositeVector*>& result)
{
  AMANZI_ASSERT(false); // not reachable, IsDifferentiableWRT() always false
}

} // namespace Relations
} // namespace Amanzi
