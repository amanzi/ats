/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Drainage rate.

*/

#include "drainage_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

DrainageEvaluator::DrainageEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  Key akey = my_keys_.front().first;
  Tag tag = my_keys_.front().second;
  Key domain = Keys::getDomain(akey);
  akey = Keys::getVarName(akey);
  my_keys_.clear();

  drainage_key_ = Keys::in(akey, "drainage") ? akey : "drainage";
  drainage_key_ = Keys::readKey(plist_, domain, "drainage", drainage_key_);
  my_keys_.emplace_back(KeyTag{ drainage_key_, tag });

  fracwet_key_ = Keys::in(akey, "fracwet") ? akey : "fracwet";
  fracwet_key_ = Keys::readKey(plist_, domain, "fraction wet", fracwet_key_);
  my_keys_.emplace_back(KeyTag{ fracwet_key_, tag });

  // Set up my dependencies.
  // -- the extent of material, LAI for
  ai_key_ = Keys::readKey(plist_, domain, "area index", "area_index");
  dependencies_.insert(KeyTag{ ai_key_, tag });

  // -- water equivalent of the layer drained, in m^3 water per m^2 grid cell area
  wc_key_ = Keys::readKey(plist_, domain, "water equivalent", "water_equivalent");
  dependencies_.insert(KeyTag{ wc_key_, tag });

  // parameters for the drainage model
  tau_ = plist_.get<double>("drainage timescale [s]", 864);

  // default from Dickinson et al 93, CLM 4.5 Tech note eqn 7.8
  // NOTE: put this in land cover!
  wc_sat_ = plist_.get<double>("saturated specific water content [m^3 H2O / m^2 leaf area]", 1.e-4);
  if (wc_sat_ < 0) {
    Errors::Message message("\"saturated specific water content\" must be greater than 0.");
    Exceptions::amanzi_throw(message);
  }
};


Teuchos::RCP<Evaluator>
DrainageEvaluator::Clone() const
{
  return Teuchos::rcp(new DrainageEvaluator(*this));
}


void
DrainageEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  Tag tag = my_keys_.front().second;
  // Pull dependencies out of state.
  const Epetra_MultiVector& wc = *S.Get<CompositeVector>(wc_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& ai = *S.Get<CompositeVector>(ai_key_, tag).ViewComponent("cell", false);

  Epetra_MultiVector& res_drainage_c = *results[0]->ViewComponent("cell", false);
  Epetra_MultiVector& res_fracwet_c = *results[1]->ViewComponent("cell", false);

  // evaluate the model
  for (int c = 0; c != res_drainage_c.MyLength(); ++c) {
    double wc_cell = std::max(wc[0][c], 0.);
    double ai_cell = std::max(ai[0][c], 0.);

    // convert from m^3 H20/ m^2 leaf area to m^3 H20 / m^2 cell area
    double wc_cell_sat = wc_sat_ * ai_cell;
    res_fracwet_c[0][c] = wc_cell_sat > 0. ? wc_cell / wc_cell_sat : 0.;

    // must be in [0,1] -- note that wc_cell can be > wc_cell_sat
    res_fracwet_c[0][c] = std::max(std::min(res_fracwet_c[0][c], 1.0), 0.0);

    if (wc_cell > wc_cell_sat) {
      //  is oversaturated and draining
      // NOTE: should this actually be:
      // res_drainage_c[0][c] = (wc_cell - wc_cell_sat) / ai[0][c] / tau_;
      // to make it proportional in units of m^3 H20 per m^2 leaf area
      res_drainage_c[0][c] = (wc_cell - wc_cell_sat) / tau_;
    } else {
      res_drainage_c[0][c] = 0.;
    }
  }
}


void
DrainageEvaluator::EvaluatePartialDerivative_(const State& S,
                                              const Key& wrt_key,
                                              const Tag& wrt_tag,
                                              const std::vector<CompositeVector*>& results)
{
  Tag tag = my_keys_.front().second;

  // Pull dependencies out of state.
  const Epetra_MultiVector& wc = *S.Get<CompositeVector>(wc_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& ai = *S.Get<CompositeVector>(ai_key_, tag).ViewComponent("cell", false);

  Epetra_MultiVector& res_drainage_c = *results[0]->ViewComponent("cell", false);
  Epetra_MultiVector& res_fracwet_c = *results[1]->ViewComponent("cell", false);

  if (wrt_key == wc_key_) {
    for (int c = 0; c != res_drainage_c.MyLength(); ++c) {
      double wc_cell = std::max(wc[0][c], 0.);
      double wc_cell_sat = wc_sat_ * ai[0][c];
      res_fracwet_c[0][c] = wc_cell_sat > 0. ? 1 / wc_cell_sat : 0;
      if (wc_cell > wc_cell_sat) {
        //  is oversaturated and draining
        res_drainage_c[0][c] = 1.0 / tau_;
      } else {
        res_drainage_c[0][c] = 0.;
      }
    }
  } else {
    AMANZI_ASSERT(0);
  }
}

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
