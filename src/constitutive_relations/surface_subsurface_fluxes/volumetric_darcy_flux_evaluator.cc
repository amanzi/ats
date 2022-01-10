#include "volumetric_darcy_flux_evaluator.hh"

namespace Amanzi {
namespace Relations {

  Volumetric_FluxEvaluator::Volumetric_FluxEvaluator(
                                                     Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {
    if (my_key_ == std::string("")) {
      my_key_ = plist_.get<std::string>("vol darcy flux key", "vol_darcy_flux");
    }

    flux_key_ = plist_.get<std::string>("flux key", "darcy_flux");
    //dependencies_.insert(flux_key_);
    dependencies_.insert("saturation_liquid");

    dens_key_ = plist_.get<std::string>("molar density key", "molar_density_liquid");
    dependencies_.insert(dens_key_);

    mesh_key_ = plist_.get<std::string>("mesh key", "domain");
  }

  Volumetric_FluxEvaluator::Volumetric_FluxEvaluator(const Volumetric_FluxEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    flux_key_(other.flux_key_),
    dens_key_(other.dens_key_),
    mesh_key_(other.mesh_key_)
  {}

  Teuchos::RCP<FieldEvaluator> Volumetric_FluxEvaluator::Clone() const {
    return Teuchos::rcp(new Volumetric_FluxEvaluator(*this));
  }


void Volumetric_FluxEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S) {
  // for now just passing... might do something later here?
  S->Require<CompositeVector,CompositeVectorSpace>(my_key_, Tags::NEXT,  my_key_);
  S->Require<CompositeVector,CompositeVectorSpace>(flux_key_, Tags::NEXT);
}

  void Volumetric_FluxEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
                                                const Teuchos::Ptr<CompositeVector>& result){

    const Epetra_MultiVector& darcy_flux = *S->Get<CompositeVector>(flux_key_).ViewComponent("face",false);
    const Epetra_MultiVector& molar_density = *S->Get<CompositeVector>(dens_key_).ViewComponent("cell",false);

    Epetra_MultiVector& res_v = *result->ViewComponent("face",false);

    Teuchos::RCP<const AmanziMesh::Mesh> mesh_ = S->GetMesh(mesh_key_);

    int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
    AmanziMesh::Entity_ID_List cells;
  
    for (int f = 0; f < nfaces_owned; f++){
      mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
      double n_liq=0.;
      for (int c=0; c<cells.size();c++) n_liq += molar_density[0][c];
      n_liq /= cells.size();
      if (n_liq > 0) res_v[0][f] = darcy_flux[0][f]/n_liq;
      else res_v[0][f] = 0.;
    }

  }


  void Volumetric_FluxEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
                                                                   Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
    AMANZI_ASSERT(0);
    // this would require differentiating flux wrt pressure, which we
    // don't do for now.
  }

}//namespace
}//namespace
