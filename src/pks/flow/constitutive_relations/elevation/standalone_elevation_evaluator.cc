/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  An elevation evaluator getting values from the volumetric mesh.

*/

#include "CompositeVectorFunctionFactory.hh"
#include "standalone_elevation_evaluator.hh"

namespace Amanzi {
namespace Flow {

StandaloneElevationEvaluator::StandaloneElevationEvaluator(Teuchos::ParameterList& plist)
  : ElevationEvaluator(plist)
{}

Teuchos::RCP<Evaluator>
StandaloneElevationEvaluator::Clone() const
{
  return Teuchos::rcp(new StandaloneElevationEvaluator(*this));
}

void
StandaloneElevationEvaluator::EvaluateElevationAndSlope_(
  const State& S,
  const std::vector<CompositeVector*>& results)
{
  Teuchos::Ptr<CompositeVector> elev = Teuchos::ptr(results[0]);
  Teuchos::Ptr<CompositeVector> slope = Teuchos::ptr(results[1]);
  Teuchos::Ptr<CompositeVector> aspect = Teuchos::ptr(results[2]);

  // If necessary, create the functions from paramater lists.
  if (elevation_function_ == Teuchos::null) {
    Teuchos::ParameterList plist = plist_.sublist("elevation function");
    std::vector<std::string> compnames;
    elevation_function_ = Functions::CreateCompositeVectorFunction(plist, elev->Map(), compnames);
    // note, should check that cells exist?
  }

  if (slope_function_ == Teuchos::null) {
    Teuchos::ParameterList slope_plist = plist_.sublist("slope function");
    std::vector<std::string> compnames;
    slope_function_ =
      Functions::CreateCompositeVectorFunction(slope_plist, slope->Map(), compnames);
    // note, should check that cells exist?
  }

  if (aspect_function_ == Teuchos::null) {
    if (plist_.isSublist("aspect function")) {
      Teuchos::ParameterList aspect_plist = plist_.sublist("aspect function");
      std::vector<std::string> compnames;
      aspect_function_ =
        Functions::CreateCompositeVectorFunction(aspect_plist, aspect->Map(), compnames);
      // note, should check that cells exist?
    }
  }

  // Evaluate the functions.
  elevation_function_->Compute(S.get_time(), elev);
  slope_function_->Compute(S.get_time(), slope);
  if (aspect_function_.get()) aspect_function_->Compute(S.get_time(), aspect);
};


} // namespace Flow
} // namespace Amanzi
