/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
ATS

Authors: Ethan Coon (ecoon@lanl.gov)

FieldEvaluator for water content.

Wrapping this conserved quantity as a field evaluator makes it easier to take
derivatives, keep updated, and the like.  The equation for this is simply:

WC = phi * (s_liquid * n_liquid + omega_gas * s_gas * n_gas)

This is simply the conserved quantity in Richards equation.
----------------------------------------------------------------------------- */


#include "interfrost_energy_evaluator.hh"

namespace Amanzi {
namespace Energy {

InterfrostEnergyEvaluator::InterfrostEnergyEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {
  my_key_ = plist_.get<std::string>("energy key", "energy");
 
  dependencies_.insert(std::string("porosity"));
  dependencies_.insert(std::string("base_porosity"));

  dependencies_.insert(std::string("saturation_liquid"));
  dependencies_.insert(std::string("molar_density_liquid"));
  dependencies_.insert(std::string("internal_energy_liquid"));
  dependencies_.insert(std::string("pressure"));

  dependencies_.insert(std::string("saturation_ice"));
  dependencies_.insert(std::string("molar_density_ice"));
  dependencies_.insert(std::string("internal_energy_ice"));

  dependencies_.insert(std::string("density_rock"));
  dependencies_.insert(std::string("internal_energy_rock"));
  //  dependencies_.insert(std::string("cell_volume"));

  //  check_derivative_ = true;

  beta_ = plist.get<double>("compressibility [1/Pa]");
};

Teuchos::RCP<FieldEvaluator>
InterfrostEnergyEvaluator::Clone() const {
  return Teuchos::rcp(new InterfrostEnergyEvaluator(*this));
};


void InterfrostEnergyEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {
  const Epetra_MultiVector& s_l = *S->Get<CompositeVector>("saturation_liquid").ViewComponent("cell",false);
  const Epetra_MultiVector& n_l = *S->Get<CompositeVector>("molar_density_liquid").ViewComponent("cell",false);
  const Epetra_MultiVector& u_l = *S->Get<CompositeVector>("internal_energy_liquid").ViewComponent("cell",false);
  const Epetra_MultiVector& pres = *S->Get<CompositeVector>("pressure").ViewComponent("cell",false);

  const Epetra_MultiVector& s_i = *S->Get<CompositeVector>("saturation_ice").ViewComponent("cell",false);
  const Epetra_MultiVector& n_i = *S->Get<CompositeVector>("molar_density_ice").ViewComponent("cell",false);
  const Epetra_MultiVector& u_i = *S->Get<CompositeVector>("internal_energy_ice").ViewComponent("cell",false);

  const Epetra_MultiVector& phi = *S->Get<CompositeVector>("porosity").ViewComponent("cell",false);
  const Epetra_MultiVector& phib = *S->Get<CompositeVector>("base_porosity").ViewComponent("cell",false);
  const Epetra_MultiVector& u_rock = *S->Get<CompositeVector>("internal_energy_rock").ViewComponent("cell",false);
  const Epetra_MultiVector& rho_rock = *S->Get<CompositeVector>("density_rock").ViewComponent("cell",false);
  const Epetra_MultiVector& cell_volume = *S->Get<CompositeVector>("cell_volume").ViewComponent("cell",false);
  Epetra_MultiVector& result_v = *result->ViewComponent("cell",false);

  int ncells = result->size("cell", false);
  for (int c=0; c!=ncells; ++c) {
    double pc = std::max(pres[0][c] - 101325., 0.);
    result_v[0][c] = phi[0][c] * (
        s_l[0][c] * n_l[0][c] * u_l[0][c] * (1+beta_*pc)
        + s_i[0][c] * n_i[0][c] * u_i[0][c])
        + (1.0 - phib[0][c]) * u_rock[0][c] * rho_rock[0][c];
    result_v[0][c] *= cell_volume[0][c];
  }
};


void InterfrostEnergyEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
  const Epetra_MultiVector& s_l = *S->Get<CompositeVector>("saturation_liquid").ViewComponent("cell",false);
  const Epetra_MultiVector& n_l = *S->Get<CompositeVector>("molar_density_liquid").ViewComponent("cell",false);
  const Epetra_MultiVector& u_l = *S->Get<CompositeVector>("internal_energy_liquid").ViewComponent("cell",false);
  const Epetra_MultiVector& pres = *S->Get<CompositeVector>("pressure").ViewComponent("cell",false);

  const Epetra_MultiVector& s_i = *S->Get<CompositeVector>("saturation_ice").ViewComponent("cell",false);
  const Epetra_MultiVector& n_i = *S->Get<CompositeVector>("molar_density_ice").ViewComponent("cell",false);
  const Epetra_MultiVector& u_i = *S->Get<CompositeVector>("internal_energy_ice").ViewComponent("cell",false);

  const Epetra_MultiVector& phi = *S->Get<CompositeVector>("porosity").ViewComponent("cell",false);
  const Epetra_MultiVector& phib = *S->Get<CompositeVector>("base_porosity").ViewComponent("cell",false);
  const Epetra_MultiVector& u_rock = *S->Get<CompositeVector>("internal_energy_rock").ViewComponent("cell",false);
  const Epetra_MultiVector& rho_rock = *S->Get<CompositeVector>("density_rock").ViewComponent("cell",false);
  const Epetra_MultiVector& cell_volume = *S->Get<CompositeVector>("cell_volume").ViewComponent("cell",false);
  Epetra_MultiVector& result_v = *result->ViewComponent("cell",false);

  if (wrt_key == "porosity") {
    for (unsigned int c=0; c!=result->size("cell"); ++c) {
      double pc = std::max(pres[0][c] - 101325., 0.);
      result_v[0][c] = s_l[0][c] * n_l[0][c] * u_l[0][c] * (1+beta_*pc)
          + s_i[0][c] * n_i[0][c] * u_i[0][c];
    }
  } else if (wrt_key == "base_porosity") {
    for (unsigned int c=0; c!=result->size("cell"); ++c) {
      result_v[0][c] = - rho_rock[0][c] * u_rock[0][c];
    }

  } else if (wrt_key == "saturation_liquid") {
    for (unsigned int c=0; c!=result->size("cell"); ++c) {
      double pc = std::max(pres[0][c] - 101325., 0.);
      result_v[0][c] = phi[0][c]*n_l[0][c]*u_l[0][c] * (1+beta_*pc);
    }
  } else if (wrt_key == "molar_density_liquid") {
    for (unsigned int c=0; c!=result->size("cell"); ++c) {
      double pc = std::max(pres[0][c] - 101325., 0.);
      result_v[0][c] = phi[0][c]*s_l[0][c]*u_l[0][c] * (1+beta_*pc);
    }
  } else if (wrt_key == "internal_energy_liquid") {
    for (unsigned int c=0; c!=result->size("cell"); ++c) {
      double pc = std::max(pres[0][c] - 101325., 0.);
      result_v[0][c] = phi[0][c]*s_l[0][c]*n_l[0][c] * (1+beta_*pc);
    }
  } else if (wrt_key == "pressure") {
    for (unsigned int c=0; c!=result->size("cell"); ++c) {
      double pc = std::max(pres[0][c] - 101325., 0.);
      result_v[0][c] = phi[0][c]*s_l[0][c]*n_l[0][c]*u_l[0][c] * beta_;
    }

  } else if (wrt_key == "saturation_ice") {
    for (unsigned int c=0; c!=result->size("cell"); ++c) {
      result_v[0][c] = phi[0][c]*n_i[0][c]*u_i[0][c];
    }
  } else if (wrt_key == "molar_density_ice") {
    for (unsigned int c=0; c!=result->size("cell"); ++c) {
      result_v[0][c] = phi[0][c]*s_i[0][c]*u_i[0][c];
    }
  } else if (wrt_key == "internal_energy_ice") {
    for (unsigned int c=0; c!=result->size("cell"); ++c) {
      result_v[0][c] = phi[0][c]*s_i[0][c]*n_i[0][c];
    }

  } else if (wrt_key == "internal_energy_rock") {
    for (unsigned int c=0; c!=result->size("cell"); ++c) {
      result_v[0][c] = (1.0 - phib[0][c])*rho_rock[0][c];
    }
  } else if (wrt_key == "density_rock") {
    for (unsigned int c=0; c!=result->size("cell"); ++c) {
      result_v[0][c] = (1.0 - phib[0][c])*u_rock[0][c];
    }
  } else {
    AMANZI_ASSERT(0);
  }

  for (unsigned int c=0; c!=result->size("cell"); ++c) {
    result_v[0][c] *= cell_volume[0][c];
  }
};


} //namespace
} //namespace
