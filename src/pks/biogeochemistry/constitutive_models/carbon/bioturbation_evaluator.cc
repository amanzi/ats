/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Bioturbation via diffusion

*/

#include "Epetra_SerialDenseVector.h"

#include "bioturbation_evaluator.hh"

namespace Amanzi {
namespace BGC {
namespace BGCRelations {

BioturbationEvaluator::BioturbationEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  Tag tag = my_keys_.front().second;
  Key domain_name = Keys::getDomain(my_keys_.front().first);

  carbon_key_ = Keys::readKey(plist_, domain_name, "soil organic matter", "soil_organic_matter");
  dependencies_.insert(KeyTag{ carbon_key_, tag });

  diffusivity_key_ =
    Keys::readKey(plist_, domain_name, "cryoturbation diffusivity", "cryoturbation_diffusivity");
  dependencies_.insert(KeyTag{ diffusivity_key_, tag });
}


Teuchos::RCP<Evaluator>
BioturbationEvaluator::Clone() const
{
  return Teuchos::rcp(new BioturbationEvaluator(*this));
}


// Required methods from EvaluatorSecondaryMonotypeCV
void
BioturbationEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  auto carbon_cv = S.GetPtr<CompositeVector>(carbon_key_, tag);
  const AmanziMesh::Mesh& mesh = *carbon_cv->Mesh();

  const Epetra_MultiVector& carbon = *carbon_cv->ViewComponent("cell", false);
  const Epetra_MultiVector& diff =
    *S.GetPtr<CompositeVector>(diffusivity_key_, tag)->ViewComponent("cell", false);
  Epetra_MultiVector& res_c = *result[0]->ViewComponent("cell", false);

  // iterate over columns of the mesh
  int ncolumns = mesh.columns.num_columns_owned;
  for (int i = 0; i < ncolumns; ++i) {
    // grab the column
    const AmanziMesh::Entity_ID_View& col = mesh.columns.getCells(i);

    Epetra_SerialDenseVector dC_up(carbon.NumVectors());
    Epetra_SerialDenseVector dC_dn(carbon.NumVectors());

    // loop over column, getting cell index ci and cell c
    int ci = 0;
    for (const auto& c :col) {
      double my_z = mesh.getCellCentroid(c)[2];
      double dz_up = 0.;
      double dz_dn = 0.;

      if (ci != 0) {
        double my_z = mesh.getCellCentroid(c)[2];
        int c_up = col[ci - 1];
        dz_up = mesh.getCellCentroid(c_up)[2] - my_z;

        for (int p = 0; p != carbon.NumVectors(); ++p) {
          dC_up[p] = (diff[p][c] + diff[p][c_up]) / 2. * (carbon[p][c] - carbon[p][c_up]) / dz_up;
        }
      }

      if (ci != col.size() - 1) {
        int c_dn = col[ci + 1];
        dz_dn = mesh.getCellCentroid(c_dn)[2] - my_z;

        for (int p = 0; p != carbon.NumVectors(); ++p) {
          dC_dn[p] = (diff[p][c] + diff[p][c_dn]) / 2. * (carbon[p][c_dn] - carbon[p][c]) / dz_up;
        }
      }

      double dz = dz_dn == 0. ? dz_up : dz_up == 0. ? dz_dn : (dz_up + dz_dn) / 2.;
      for (int p = 0; p != carbon.NumVectors(); ++p) { res_c[p][c] = (dC_dn[p] - dC_up[p]) / dz; }
      ++ci; 
    }
  }
}


void
BioturbationEvaluator::EvaluatePartialDerivative_(const State& S,
                                                  const Key& wrt_key,
                                                  const Tag& wrt_tag,
                                                  const std::vector<CompositeVector*>& result)
{
  AMANZI_ASSERT(0);
}


} // namespace BGCRelations
} // namespace BGC
} // namespace Amanzi
