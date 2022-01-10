/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
ATS

Authors: Ethan Coon (ecoon@lanl.gov)

FieldEvaluator for water content.

INTERFROST's comparison uses a very odd compressibility term that doesn't
quite fit into either compressible porosity or into a compressible density, so
it needs a special evaluator.
----------------------------------------------------------------------------- */


#include "interfrost_water_content.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

InterfrostWaterContent::InterfrostWaterContent(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {
  my_key_ = std::string("water_content");

  dependencies_.insert(std::string("porosity"));

  dependencies_.insert(std::string("saturation_liquid"));
  dependencies_.insert(std::string("molar_density_liquid"));

  dependencies_.insert(std::string("saturation_ice"));
  dependencies_.insert(std::string("molar_density_ice"));

  dependencies_.insert(std::string("pressure"));

  beta_ = plist.get<double>("compressibility [1/Pa]");
  
  //  dependencies_.insert(std::string("cell_volume"));

  //  check_derivative_ = true;
};

Teuchos::RCP<FieldEvaluator>
InterfrostWaterContent::Clone() const {
  return Teuchos::rcp(new InterfrostWaterContent(*this));
};


void InterfrostWaterContent::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {
  const Epetra_MultiVector& s_l = *S->Get<CompositeVector>("saturation_liquid").ViewComponent("cell",false);
  const Epetra_MultiVector& n_l = *S->Get<CompositeVector>("molar_density_liquid").ViewComponent("cell",false);

  const Epetra_MultiVector& s_i = *S->Get<CompositeVector>("saturation_ice").ViewComponent("cell",false);
  const Epetra_MultiVector& n_i = *S->Get<CompositeVector>("molar_density_ice").ViewComponent("cell",false);

  const Epetra_MultiVector& phi = *S->Get<CompositeVector>("porosity").ViewComponent("cell",false);
  const Epetra_MultiVector& pressure = *S->Get<CompositeVector>("pressure").ViewComponent("cell",false);
  const Epetra_MultiVector& cell_volume = *S->Get<CompositeVector>("cell_volume").ViewComponent("cell",false);
  Epetra_MultiVector& result_v = *result->ViewComponent("cell",false);

  int ncells = result->size("cell",false);
  for (int c=0; c!=ncells; ++c) {
    double pr = std::max(pressure[0][c] - 101325., 0.);
    result_v[0][c] = phi[0][c] * ( s_l[0][c]*n_l[0][c]*(1+beta_*pr)
            + s_i[0][c]*n_i[0][c]);
    result_v[0][c] *= cell_volume[0][c];
  }
};


void InterfrostWaterContent::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
  const Epetra_MultiVector& s_l = *S->Get<CompositeVector>("saturation_liquid").ViewComponent("cell",false);
  const Epetra_MultiVector& n_l = *S->Get<CompositeVector>("molar_density_liquid").ViewComponent("cell",false);

  const Epetra_MultiVector& s_i = *S->Get<CompositeVector>("saturation_ice").ViewComponent("cell",false);
  const Epetra_MultiVector& n_i = *S->Get<CompositeVector>("molar_density_ice").ViewComponent("cell",false);

  const Epetra_MultiVector& pressure = *S->Get<CompositeVector>("pressure").ViewComponent("cell",false);
  const Epetra_MultiVector& phi = *S->Get<CompositeVector>("porosity").ViewComponent("cell",false);
  const Epetra_MultiVector& cell_volume = *S->Get<CompositeVector>("cell_volume").ViewComponent("cell",false);
  Epetra_MultiVector& result_v = *result->ViewComponent("cell",false);

  int ncells = result->size("cell",false);
  if (wrt_key == "porosity") {
    for (int c=0; c!=ncells; ++c) {
      double pr = std::max(pressure[0][c] - 101325., 0.);
      result_v[0][c] = ( s_l[0][c]*n_l[0][c]*(1+beta_*pr)
              + s_i[0][c]*n_i[0][c]);
    }
  } else if (wrt_key == "saturation_liquid") {
    for (int c=0; c!=ncells; ++c) {
      double pr = std::max(pressure[0][c] - 101325., 0.);
      result_v[0][c] = phi[0][c] * n_l[0][c] * (1+beta_*pr);
    }
  } else if (wrt_key == "molar_density_liquid") {
    for (int c=0; c!=ncells; ++c) {
      double pr = std::max(pressure[0][c] - 101325., 0.);
      result_v[0][c] = phi[0][c] * s_l[0][c] * (1+beta_*pr);
    }
  } else if (wrt_key == "saturation_ice") {
    for (int c=0; c!=ncells; ++c) {
      result_v[0][c] = phi[0][c] * n_i[0][c];
    }
  } else if (wrt_key == "molar_density_ice") {
    for (int c=0; c!=ncells; ++c) {
      result_v[0][c] = phi[0][c] * s_i[0][c];
    }
  } else if (wrt_key == "pressure") {
    for (int c=0; c!=ncells; ++c) {
      result_v[0][c] = phi[0][c] * s_l[0][c] * n_l[0][c] * beta_;
    }
  }

  for (int c=0; c!=ncells; ++c) {
    result_v[0][c] *= cell_volume[0][c];
  }
};


} //namespace
} //namespace
} //namespace
