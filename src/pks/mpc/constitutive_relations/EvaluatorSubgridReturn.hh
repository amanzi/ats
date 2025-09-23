/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/
/*!

Computes the subgrid return from the hyporheic zone back to the surface system,
which acts like a source term.


.. math::
   Q = \alpha \Theta \int_0^\infty \xi d\tau

Note that the :math:`\tau` are equally spaced in travel time, so this integral
simplifies to a summation over grid cells in the subgrid domain, divided by the
number of grid cells.

`"evaluator type`" = `"subgrid return`"

.. _evaluator-subgrid-return-spec:
.. admonition:: evaluator-subgrid-return-spec

   KEYS:
   - `"exchange coefficient`"
   - `"liquid water content`" **water_content**
   - `"exchanged key`" **mole_fraction**

*/

#pragma once

#include <string>

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace ATS_Physics {

class EvaluatorSubgridReturn : public EvaluatorSecondaryMonotypeCV {
 public:
  EvaluatorSubgridReturn(Teuchos::ParameterList& plist);
  EvaluatorSubgridReturn(const EvaluatorSubgridReturn& other) = default;

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override
  {
    AMANZI_ASSERT(false);
  }

  virtual bool IsDifferentiableWRT(const State& S, const Key& key, const Tag& tag) const override
  {
    return false;
  }

  virtual void EnsureCompatibility_ToDeps_(State& S) override;

 private:
  int n_dofs_;
  Key alpha_key_;
  Key mf_suffix_;
  Key lwc_key_;
  Key cv_key_;

  std::string domain_;
  std::string domain_set_;

  static Utils::RegisteredFactory<Evaluator, EvaluatorSubgridReturn> fac_;
};

} // namespace ATS_Physics
} // namespace Amanzi
