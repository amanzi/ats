/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Interface for a thermal conductivity model with two phases.

*/

#include "thermal_conductivity_surface_evaluator.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Energy {

ThermalConductivitySurfaceEvaluator::ThermalConductivitySurfaceEvaluator(
  Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  Key domain = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;
  uf_key_ = Keys::readKey(plist_, domain, "unfrozen fraction", "unfrozen_fraction");
  dependencies_.insert(KeyTag{ uf_key_, tag });

  height_key_ = Keys::readKey(plist_, domain, "ponded depth", "ponded_depth");
  dependencies_.insert(KeyTag{ height_key_, tag });

  AMANZI_ASSERT(plist_.isSublist("thermal conductivity parameters"));
  Teuchos::ParameterList sublist = plist_.sublist("thermal conductivity parameters");
  K_liq_ = sublist.get<double>("thermal conductivity of water [W m^-1 K^-1]", 0.58);
  K_ice_ = sublist.get<double>("thermal conductivity of ice [W m^-1 K^-1]", 2.18);
  min_K_ = sublist.get<double>("minimum thermal conductivity", 1.e-14);
}

Teuchos::RCP<Evaluator>
ThermalConductivitySurfaceEvaluator::Clone() const
{
  return Teuchos::rcp(new ThermalConductivitySurfaceEvaluator(*this));
}

void
ThermalConductivitySurfaceEvaluator::Evaluate_(const State& S,
                                               const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  // pull out the dependencies
  Teuchos::RCP<const CompositeVector> eta = S.GetPtr<CompositeVector>(uf_key_, tag);
  Teuchos::RCP<const CompositeVector> height = S.GetPtr<CompositeVector>(height_key_, tag);

  for (CompositeVector::name_iterator comp = result[0]->begin(); comp != result[0]->end(); ++comp) {
    // much more efficient to pull out vectors first
    const Epetra_MultiVector& eta_v = *eta->ViewComponent(*comp, false);
    const Epetra_MultiVector& height_v = *height->ViewComponent(*comp, false);
    Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

    int ncomp = result[0]->size(*comp, false);
    for (int i = 0; i != ncomp; ++i) {
      result_v[0][i] =
        std::max(min_K_, height_v[0][i] * (K_liq_ * eta_v[0][i] + K_ice_ * (1. - eta_v[0][i])));
    }
  }

  result[0]->Scale(1.e-6); // convert to MJ
}


void
ThermalConductivitySurfaceEvaluator::EvaluatePartialDerivative_(
  const State& S,
  const Key& wrt_key,
  const Tag& wrt_tag,
  const std::vector<CompositeVector*>& result)
{
  std::cout << "THERMAL CONDUCITIVITY: Derivative not implemented yet!" << wrt_key << "\n";
  AMANZI_ASSERT(0);        // not implemented, not yet needed
  result[0]->Scale(1.e-6); // convert to MJ
}

} // namespace Energy
} // namespace ATS_Physics
} // namespace Amanzi
