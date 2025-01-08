/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/* -----------------------------------------------------------------------------
ATS

Evaluator for water content.

Wrapping this conserved quantity as a field evaluator makes it easier to take
derivatives, keep updated, and the like.  The equation for this is simply:

WC = phi * (s_liquid * n_liquid + omega_gas * s_gas * n_gas)

This is simply the conserved quantity in Richards equation.
----------------------------------------------------------------------------- */


#include "interfrost_energy_evaluator.hh"

namespace Amanzi {
namespace Energy {

InterfrostEnergyEvaluator::InterfrostEnergyEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  Key domain_name = Keys::getDomain(my_keys_.front().first);
  auto tag = my_keys_.front().second;

  dependencies_.insert(KeyTag{ std::string("porosity"), tag });
  dependencies_.insert(KeyTag{ std::string("base_porosity"), tag });

  dependencies_.insert(KeyTag{ std::string("saturation_liquid"), tag });
  dependencies_.insert(KeyTag{ std::string("molar_density_liquid"), tag });
  dependencies_.insert(KeyTag{ std::string("internal_energy_liquid"), tag });
  dependencies_.insert(KeyTag{ std::string("pressure"), tag });

  dependencies_.insert(KeyTag{ std::string("saturation_ice"), tag });
  dependencies_.insert(KeyTag{ std::string("molar_density_ice"), tag });
  dependencies_.insert(KeyTag{ std::string("internal_energy_ice"), tag });

  dependencies_.insert(KeyTag{ std::string("density_rock"), tag });
  dependencies_.insert(KeyTag{ std::string("internal_energy_rock"), tag });
  //  dependencies_.insert(std::string("cell_volume"));

  //  check_derivative_ = true;

  beta_ = plist.get<double>("compressibility [1/Pa]");
};

Teuchos::RCP<Evaluator>
InterfrostEnergyEvaluator::Clone() const
{
  return Teuchos::rcp(new InterfrostEnergyEvaluator(*this));
};


void
InterfrostEnergyEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  auto tag = my_keys_.front().second;
  const Epetra_MultiVector& s_l =
    *S.Get<CompositeVector>("saturation_liquid", tag).ViewComponent("cell", false);
  const Epetra_MultiVector& n_l =
    *S.Get<CompositeVector>("molar_density_liquid", tag).ViewComponent("cell", false);
  const Epetra_MultiVector& u_l =
    *S.Get<CompositeVector>("internal_energy_liquid", tag).ViewComponent("cell", false);
  const Epetra_MultiVector& pres =
    *S.Get<CompositeVector>("pressure", tag).ViewComponent("cell", false);

  const Epetra_MultiVector& s_i =
    *S.Get<CompositeVector>("saturation_ice", tag).ViewComponent("cell", false);
  const Epetra_MultiVector& n_i =
    *S.Get<CompositeVector>("molar_density_ice", tag).ViewComponent("cell", false);
  const Epetra_MultiVector& u_i =
    *S.Get<CompositeVector>("internal_energy_ice", tag).ViewComponent("cell", false);

  const Epetra_MultiVector& phi =
    *S.Get<CompositeVector>("porosity", tag).ViewComponent("cell", false);
  const Epetra_MultiVector& phib =
    *S.Get<CompositeVector>("base_porosity", tag).ViewComponent("cell", false);
  const Epetra_MultiVector& u_rock =
    *S.Get<CompositeVector>("internal_energy_rock", tag).ViewComponent("cell", false);
  const Epetra_MultiVector& rho_rock =
    *S.Get<CompositeVector>("density_rock", tag).ViewComponent("cell", false);
  const Epetra_MultiVector& cell_volume =
    *S.Get<CompositeVector>("cell_volume", tag).ViewComponent("cell", false);
  Epetra_MultiVector& result_v = *result[0]->ViewComponent("cell", false);

  int ncells = result[0]->size("cell", false);
  for (int c = 0; c != ncells; ++c) {
    double pc = std::max(pres[0][c] - 101325., 0.);
    result_v[0][c] = phi[0][c] * (s_l[0][c] * n_l[0][c] * u_l[0][c] * (1 + beta_ * pc) +
                                  s_i[0][c] * n_i[0][c] * u_i[0][c]) +
                     (1.0 - phib[0][c]) * u_rock[0][c] * rho_rock[0][c];
    result_v[0][c] *= cell_volume[0][c];
  }
};


void
InterfrostEnergyEvaluator::EvaluatePartialDerivative_(const State& S,
                                                      const Key& wrt_key,
                                                      const Tag& wrt_tag,
                                                      const std::vector<CompositeVector*>& result)
{
  auto tag = my_keys_.front().second;
  const Epetra_MultiVector& s_l =
    *S.Get<CompositeVector>("saturation_liquid", tag).ViewComponent("cell", false);
  const Epetra_MultiVector& n_l =
    *S.Get<CompositeVector>("molar_density_liquid", tag).ViewComponent("cell", false);
  const Epetra_MultiVector& u_l =
    *S.Get<CompositeVector>("internal_energy_liquid", tag).ViewComponent("cell", false);
  const Epetra_MultiVector& pres =
    *S.Get<CompositeVector>("pressure", tag).ViewComponent("cell", false);

  const Epetra_MultiVector& s_i =
    *S.Get<CompositeVector>("saturation_ice", tag).ViewComponent("cell", false);
  const Epetra_MultiVector& n_i =
    *S.Get<CompositeVector>("molar_density_ice", tag).ViewComponent("cell", false);
  const Epetra_MultiVector& u_i =
    *S.Get<CompositeVector>("internal_energy_ice", tag).ViewComponent("cell", false);

  const Epetra_MultiVector& phi =
    *S.Get<CompositeVector>("porosity", tag).ViewComponent("cell", false);
  const Epetra_MultiVector& phib =
    *S.Get<CompositeVector>("base_porosity", tag).ViewComponent("cell", false);
  const Epetra_MultiVector& u_rock =
    *S.Get<CompositeVector>("internal_energy_rock", tag).ViewComponent("cell", false);
  const Epetra_MultiVector& rho_rock =
    *S.Get<CompositeVector>("density_rock", tag).ViewComponent("cell", false);
  const Epetra_MultiVector& cell_volume =
    *S.Get<CompositeVector>("cell_volume", tag).ViewComponent("cell", false);
  Epetra_MultiVector& result_v = *result[0]->ViewComponent("cell", false);

  if (wrt_key == "porosity") {
    for (unsigned int c = 0; c != result[0]->size("cell"); ++c) {
      double pc = std::max(pres[0][c] - 101325., 0.);
      result_v[0][c] =
        s_l[0][c] * n_l[0][c] * u_l[0][c] * (1 + beta_ * pc) + s_i[0][c] * n_i[0][c] * u_i[0][c];
    }
  } else if (wrt_key == "base_porosity") {
    for (unsigned int c = 0; c != result[0]->size("cell"); ++c) {
      result_v[0][c] = -rho_rock[0][c] * u_rock[0][c];
    }

  } else if (wrt_key == "saturation_liquid") {
    for (unsigned int c = 0; c != result[0]->size("cell"); ++c) {
      double pc = std::max(pres[0][c] - 101325., 0.);
      result_v[0][c] = phi[0][c] * n_l[0][c] * u_l[0][c] * (1 + beta_ * pc);
    }
  } else if (wrt_key == "molar_density_liquid") {
    for (unsigned int c = 0; c != result[0]->size("cell"); ++c) {
      double pc = std::max(pres[0][c] - 101325., 0.);
      result_v[0][c] = phi[0][c] * s_l[0][c] * u_l[0][c] * (1 + beta_ * pc);
    }
  } else if (wrt_key == "internal_energy_liquid") {
    for (unsigned int c = 0; c != result[0]->size("cell"); ++c) {
      double pc = std::max(pres[0][c] - 101325., 0.);
      result_v[0][c] = phi[0][c] * s_l[0][c] * n_l[0][c] * (1 + beta_ * pc);
    }
  } else if (wrt_key == "pressure") {
    for (unsigned int c = 0; c != result[0]->size("cell"); ++c) {
      double pc = std::max(pres[0][c] - 101325., 0.);
      result_v[0][c] = phi[0][c] * s_l[0][c] * n_l[0][c] * u_l[0][c] * beta_;
    }

  } else if (wrt_key == "saturation_ice") {
    for (unsigned int c = 0; c != result[0]->size("cell"); ++c) {
      result_v[0][c] = phi[0][c] * n_i[0][c] * u_i[0][c];
    }
  } else if (wrt_key == "molar_density_ice") {
    for (unsigned int c = 0; c != result[0]->size("cell"); ++c) {
      result_v[0][c] = phi[0][c] * s_i[0][c] * u_i[0][c];
    }
  } else if (wrt_key == "internal_energy_ice") {
    for (unsigned int c = 0; c != result[0]->size("cell"); ++c) {
      result_v[0][c] = phi[0][c] * s_i[0][c] * n_i[0][c];
    }

  } else if (wrt_key == "internal_energy_rock") {
    for (unsigned int c = 0; c != result[0]->size("cell"); ++c) {
      result_v[0][c] = (1.0 - phib[0][c]) * rho_rock[0][c];
    }
  } else if (wrt_key == "density_rock") {
    for (unsigned int c = 0; c != result[0]->size("cell"); ++c) {
      result_v[0][c] = (1.0 - phib[0][c]) * u_rock[0][c];
    }
  } else {
    AMANZI_ASSERT(0);
  }

  for (unsigned int c = 0; c != result[0]->size("cell"); ++c) {
    result_v[0][c] *= cell_volume[0][c];
  }
};


} // namespace Energy
} // namespace Amanzi
