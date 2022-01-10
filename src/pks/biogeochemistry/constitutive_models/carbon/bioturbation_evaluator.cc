/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Bioturbation via diffusion

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "Epetra_SerialDenseVector.h"

#include "bioturbation_evaluator.hh"

namespace Amanzi {
namespace BGC {
namespace BGCRelations {

BioturbationEvaluator::BioturbationEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {

  carbon_key_ = plist_.get<std::string>("SOM key", "soil_organic_matter");
  dependencies_.insert(carbon_key_);
  diffusivity_key_ = plist_.get<std::string>("cryoturbation diffusivity key", "cryoturbation_diffusivity");
  dependencies_.insert(diffusivity_key_);

  if (my_key_ == std::string("")) {
    my_key_ = plist_.get<std::string>("divergence of bioturbation fluxes",
            "div_bioturbation");
  }
}


BioturbationEvaluator::BioturbationEvaluator(const BioturbationEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    carbon_key_(other.carbon_key_),
    diffusivity_key_(other.diffusivity_key_) {}

Teuchos::RCP<FieldEvaluator>
BioturbationEvaluator::Clone() const {
  return Teuchos::rcp(new BioturbationEvaluator(*this));
}


// Required methods from SecondaryVariableFieldEvaluator
void BioturbationEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {

  Teuchos::RCP<const CompositeVector> carbon_cv = S->GetPtr<CompositeVector>(carbon_key_);
  const AmanziMesh::Mesh& mesh = *carbon_cv->Mesh();
  
  const Epetra_MultiVector& carbon = *carbon_cv->ViewComponent("cell",false);
  const Epetra_MultiVector& diff = *S->GetPtr<CompositeVector>(diffusivity_key_)
      ->ViewComponent("cell",false);
  Epetra_MultiVector& res_c = *result->ViewComponent("cell",false);

  // iterate over columns of the mesh
  int ncolumns = mesh.num_columns();
  for (int i=0; i<ncolumns; ++i) {
    // grab the column
    const AmanziMesh::Entity_ID_List& col = mesh.cells_of_column(i);

    Epetra_SerialDenseVector dC_up(carbon.NumVectors());
    Epetra_SerialDenseVector dC_dn(carbon.NumVectors());

    // loop over column, getting cell index ci and cell c
    int ci=0;
    for (AmanziMesh::Entity_ID_List::const_iterator c=col.begin(); c!=col.end(); ++c, ++ci) {
      double my_z = mesh.cell_centroid(*c)[2];
      double dz_up = 0.;
      double dz_dn = 0.;
      
      if (ci != 0) {
        double my_z = mesh.cell_centroid(*c)[2];
        int c_up = col[ci-1];
        dz_up = mesh.cell_centroid(c_up)[2] - my_z;

        for (int p=0; p!=carbon.NumVectors(); ++p) {
          dC_up[p] = (diff[p][*c]+diff[p][c_up]) / 2. * (carbon[p][*c] - carbon[p][c_up]) / dz_up;
        }
      }

      if (ci != col.size()-1) {
        int c_dn = col[ci+1];
        dz_dn = mesh.cell_centroid(c_dn)[2] - my_z;

        for (int p=0; p!=carbon.NumVectors(); ++p) {
          dC_dn[p] = (diff[p][*c]+diff[p][c_dn]) / 2. * (carbon[p][c_dn] - carbon[p][*c]) / dz_up;
        }
      }

      double dz = dz_dn == 0. ? dz_up :
          dz_up == 0. ? dz_dn : (dz_up + dz_dn) / 2.;
      for (int p=0; p!=carbon.NumVectors(); ++p) {
        res_c[p][*c] = (dC_dn[p] - dC_up[p]) / dz;
      }
    }
    
  }
}


void BioturbationEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
  AMANZI_ASSERT(0);
}



} //namespace
} //namespace
} //namespace
