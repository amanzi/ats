/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Downregulates bare soil evaporation through a dessicated zone via soil resistance.
*/


#include "evaporation_downregulation_evaluator.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace SurfaceBalance {
namespace Relations {

// Constructor from ParameterList
EvaporationDownregulationEvaluator::EvaporationDownregulationEvaluator(
  Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  Tag tag = my_keys_.front().second;
  Key domain_surf_ = Keys::getDomain(my_keys_.front().first);

  // dependency: potential_evaporation and soil resistance on surface
  pot_evap_key_ =
    Keys::readKey(plist_, domain_surf_, "potential evaporation", "potential_evaporation");
  dependencies_.insert(KeyTag{ pot_evap_key_, tag });
  rsoil_key_ = Keys::readKey(plist_, domain_surf_, "soil resistance", "soil_resistance");
  dependencies_.insert(KeyTag{ rsoil_key_, tag });
}

// Virtual copy constructor
Teuchos::RCP<Evaluator>
EvaporationDownregulationEvaluator::Clone() const
{
  return Teuchos::rcp(new EvaporationDownregulationEvaluator(*this));
}


void
EvaporationDownregulationEvaluator::Evaluate_(const State& S,
                                              const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  const Epetra_MultiVector& rsoil =
    *S.Get<CompositeVector>(rsoil_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& pot_evap =
    *S.Get<CompositeVector>(pot_evap_key_, tag).ViewComponent("cell", false);
  Epetra_MultiVector& surf_evap = *result[0]->ViewComponent("cell", false);

  for (int sc = 0; sc != surf_evap.MyLength(); ++sc) {
    surf_evap[0][sc] = pot_evap[0][sc] / (1. + rsoil[0][sc]);
  }
}


void
EvaporationDownregulationEvaluator::EvaluatePartialDerivative_(
  const State& S,
  const Key& wrt_key,
  const Tag& wrt_tag,
  const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  if (wrt_key == pot_evap_key_) {
    const Epetra_MultiVector& rsoil =
      *S.Get<CompositeVector>(rsoil_key_, tag).ViewComponent("cell", false);
    const Epetra_MultiVector& pot_evap =
      *S.Get<CompositeVector>(pot_evap_key_, tag).ViewComponent("cell", false);
    Epetra_MultiVector& surf_evap = *result[0]->ViewComponent("cell", false);

    for (int sc = 0; sc != surf_evap.MyLength(); ++sc) {
      surf_evap[0][sc] = 1. / (1. + rsoil[0][sc]);
    }
  } else if (wrt_key == rsoil_key_) {
    const Epetra_MultiVector& rsoil =
      *S.Get<CompositeVector>(rsoil_key_, tag).ViewComponent("cell", false);
    const Epetra_MultiVector& pot_evap =
      *S.Get<CompositeVector>(pot_evap_key_, tag).ViewComponent("cell", false);
    Epetra_MultiVector& surf_evap = *result[0]->ViewComponent("cell", false);

    for (int sc = 0; sc != surf_evap.MyLength(); ++sc) {
      surf_evap[0][sc] = -pot_evap[0][sc] * std::pow(1. + rsoil[0][sc], -2);
    }
  } else {
    AMANZI_ASSERT(0);
  }
}


} // namespace Relations
} // namespace SurfaceBalance
} // namespace ATS_Physics
} // namespace Amanzi
