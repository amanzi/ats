/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

//! An evaluator for rerouting precip due to impervious surface.

/*!

This evaluator reroutes incoming precipitation, taking a portion of it (where
the portion is determined by the impervious area fraction) and moving it
(instantly) into the nearby stream network.

Note: this assumes that the runoff reciever is constant in time!

.. _impervious-interception-evaluator-spec:
.. admonition:: impervious-interception-evaluator-spec

   * `"maximum specific diversion rate [m s^-1]`" ``[double]`` **inf**
     Maximum rate of water removal through storm drains, etc, in units of m^3
     water per second per m^2 of _impervious_ area (specific area).

   KEYS:
   - `"impervious fraction`" **DOMAIN-impervious_fraction** The fraction of
     surface area that is impervious, this also defines the fraction of precip
     that is rerouted.
   - `"impervious runoff receiver`" **DOMAIN-impervious_runoff_receiver`" The
     Global ID of the cell that will recieve water from this cell.
   - `"incoming water source`" **DOMAIN-precipitation_rain** The source of
     water to be re-reouted -- this is typically rain, but might be
     canopy-throughfall_drainage_rain, and might be snow-melt, etc.
   - `"cell volume`" **DOMAIN-cell_volume**

*/

#pragma once

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "FunctionFactory.hh"

class Epetra_Import;

namespace Amanzi {
namespace Flow {
namespace Relations {

class ImperviousInterceptionEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit ImperviousInterceptionEvaluator(Teuchos::ParameterList& plist);
  ImperviousInterceptionEvaluator(const ImperviousInterceptionEvaluator& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new ImperviousInterceptionEvaluator(*this));
  }

  virtual bool
  IsDifferentiableWRT(const State& S, const Key& wrt_key, const Tag& wrt_tag) const override
  {
    // calculate of derivatives of this is a tricky thing to do, with
    // non-cell-local terms due to rescaling.  Just turn off derivatives
    // instead.
    return false;
  }

 protected:

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                        const Key& wrt_key,
                                        const Tag& wrt_tag,
                                        const std::vector<CompositeVector*>& result) override {};

 protected:
  Key imp_frac_key_;
  Key imp_rec_id_key_;
  Key src_key_;
  Key cv_key_;
  double Qs_max_;

  Teuchos::RCP<Epetra_Import> importer_;

 private:
  static Utils::RegisteredFactory<Evaluator, ImperviousInterceptionEvaluator> reg_;
};

} // namespace Relations
} // namespace Flow
} // namespace Amanzi

