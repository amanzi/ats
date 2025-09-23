/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/* -----------------------------------------------------------------------------
ATS

Evaluator for water content.

INTERFROST's comparison uses a very odd compressibility term that doesn't
quite fit into either compressible porosity or into a compressible density, so
it needs a special evaluator.
----------------------------------------------------------------------------- */


#include "interfrost_water_content.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {
namespace Relations {

InterfrostWaterContent::InterfrostWaterContent(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  Tag tag = my_keys_.front().second;

  dependencies_.insert(KeyTag{ std::string("porosity"), tag });
  dependencies_.insert(KeyTag{ std::string("saturation_liquid"), tag });
  dependencies_.insert(KeyTag{ std::string("molar_density_liquid"), tag });
  dependencies_.insert(KeyTag{ std::string("saturation_ice"), tag });
  dependencies_.insert(KeyTag{ std::string("molar_density_ice"), tag });
  dependencies_.insert(KeyTag{ std::string("pressure"), tag });

  beta_ = plist.get<double>("compressibility [1/Pa]");
};

Teuchos::RCP<Evaluator>
InterfrostWaterContent::Clone() const
{
  return Teuchos::rcp(new InterfrostWaterContent(*this));
};


void
InterfrostWaterContent::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  const Epetra_MultiVector& s_l =
    *S.Get<CompositeVector>("saturation_liquid", tag).ViewComponent("cell", false);
  const Epetra_MultiVector& n_l =
    *S.Get<CompositeVector>("molar_density_liquid", tag).ViewComponent("cell", false);

  const Epetra_MultiVector& s_i =
    *S.Get<CompositeVector>("saturation_ice", tag).ViewComponent("cell", false);
  const Epetra_MultiVector& n_i =
    *S.Get<CompositeVector>("molar_density_ice", tag).ViewComponent("cell", false);

  const Epetra_MultiVector& phi =
    *S.Get<CompositeVector>("porosity", tag).ViewComponent("cell", false);
  const Epetra_MultiVector& pressure =
    *S.Get<CompositeVector>("pressure", tag).ViewComponent("cell", false);
  const Epetra_MultiVector& cell_volume =
    *S.Get<CompositeVector>("cell_volume", tag).ViewComponent("cell", false);
  Epetra_MultiVector& result_v = *result[0]->ViewComponent("cell", false);

  int ncells = result[0]->size("cell", false);
  for (int c = 0; c != ncells; ++c) {
    double pr = std::max(pressure[0][c] - 101325., 0.);
    result_v[0][c] = phi[0][c] * (s_l[0][c] * n_l[0][c] * (1 + beta_ * pr) + s_i[0][c] * n_i[0][c]);
    result_v[0][c] *= cell_volume[0][c];
  }
};


void
InterfrostWaterContent::EvaluatePartialDerivative_(const State& S,
                                                   const Key& wrt_key,
                                                   const Tag& wrt_tag,
                                                   const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  const Epetra_MultiVector& s_l =
    *S.Get<CompositeVector>("saturation_liquid", tag).ViewComponent("cell", false);
  const Epetra_MultiVector& n_l =
    *S.Get<CompositeVector>("molar_density_liquid", tag).ViewComponent("cell", false);

  const Epetra_MultiVector& s_i =
    *S.Get<CompositeVector>("saturation_ice", tag).ViewComponent("cell", false);
  const Epetra_MultiVector& n_i =
    *S.Get<CompositeVector>("molar_density_ice", tag).ViewComponent("cell", false);

  const Epetra_MultiVector& pressure =
    *S.Get<CompositeVector>("pressure", tag).ViewComponent("cell", false);
  const Epetra_MultiVector& phi =
    *S.Get<CompositeVector>("porosity", tag).ViewComponent("cell", false);
  const Epetra_MultiVector& cell_volume =
    *S.Get<CompositeVector>("cell_volume", tag).ViewComponent("cell", false);
  Epetra_MultiVector& result_v = *result[0]->ViewComponent("cell", false);

  int ncells = result[0]->size("cell", false);
  if (wrt_key == "porosity") {
    for (int c = 0; c != ncells; ++c) {
      double pr = std::max(pressure[0][c] - 101325., 0.);
      result_v[0][c] = (s_l[0][c] * n_l[0][c] * (1 + beta_ * pr) + s_i[0][c] * n_i[0][c]);
    }
  } else if (wrt_key == "saturation_liquid") {
    for (int c = 0; c != ncells; ++c) {
      double pr = std::max(pressure[0][c] - 101325., 0.);
      result_v[0][c] = phi[0][c] * n_l[0][c] * (1 + beta_ * pr);
    }
  } else if (wrt_key == "molar_density_liquid") {
    for (int c = 0; c != ncells; ++c) {
      double pr = std::max(pressure[0][c] - 101325., 0.);
      result_v[0][c] = phi[0][c] * s_l[0][c] * (1 + beta_ * pr);
    }
  } else if (wrt_key == "saturation_ice") {
    for (int c = 0; c != ncells; ++c) {
      result_v[0][c] = phi[0][c] * n_i[0][c];
    }
  } else if (wrt_key == "molar_density_ice") {
    for (int c = 0; c != ncells; ++c) {
      result_v[0][c] = phi[0][c] * s_i[0][c];
    }
  } else if (wrt_key == "pressure") {
    for (int c = 0; c != ncells; ++c) {
      result_v[0][c] = phi[0][c] * s_l[0][c] * n_l[0][c] * beta_;
    }
  }

  for (int c = 0; c != ncells; ++c) {
    result_v[0][c] *= cell_volume[0][c];
  }
};


} // namespace Relations
} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi
