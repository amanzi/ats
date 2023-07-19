/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Transport PK

*/

#include "Teuchos_RCP.hpp"
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"

#include "MFD3D_Diffusion.hh"
#include "nlfv.hh"
#include "Tensor.hh"

#include "TransportDefs.hh"
#include "transport_ats.hh"

namespace Amanzi {
namespace Transport {

/* *******************************************************************
* Calculate dispersive tensor from given water fluxes. The flux is
* assumed to be scaled by face area.
******************************************************************* */
void
Transport_ATS::CalculateDispersionTensor_(const Epetra_MultiVector& water_flux,
                                          const Epetra_MultiVector& porosity,
                                          const Epetra_MultiVector& saturation,
                                          const Epetra_MultiVector& mol_density)
{
  D_.resize(ncells_owned);
  for (int c = 0; c < ncells_owned; c++) D_[c].Init(dim, 1);

  AmanziGeometry::Point velocity(dim);
  WhetStone::MFD3D_Diffusion mfd3d(mesh_);
  WhetStone::Polynomial poly(dim, 1);

  for (int c = 0; c < ncells_owned; ++c) {
    auto faces = mesh_->getCellFaces(c);
    int nfaces = faces.size();

    std::vector<WhetStone::Polynomial> flux(nfaces);
    for (int n = 0; n < nfaces; n++) {
      flux[n].Reshape(dim, 0);
      flux[n](0) = water_flux[0][faces[n]];
    }
    mfd3d.L2Cell(c, flux, flux, NULL, poly);

    // note that while this is called velocity, it is actually a vector-valued,
    // cell-centered water mass flux.  There is no division by density.  That
    // means the factor of density should not be required here? --ETC
    for (int k = 0; k < dim; ++k) velocity[k] = poly(k + 1);
    D_[c] = mdm_->second[(*mdm_->first)[c]]->mech_dispersion(
      velocity, axi_symmetry_[c], saturation[0][c], porosity[0][c]);
    // double mol_den = mol_density[0][c];
    // D_[c] *= mol_den;
  }
}


/* *******************************************************************
* Calculate diffusion tensor and add it to the dispersion tensor.
******************************************************************* */
void
Transport_ATS::CalculateDiffusionTensor_(double md,
                                         int phase,
                                         const Epetra_MultiVector& porosity,
                                         const Epetra_MultiVector& saturation,
                                         const Epetra_MultiVector& mol_density)
{
  if (D_.size() == 0) {
    D_.resize(ncells_owned);
    for (int c = 0; c < ncells_owned; c++) D_[c].Init(dim, 1);
  }

  for (int mb = 0; mb < mat_properties_.size(); mb++) {
    Teuchos::RCP<MaterialProperties> spec = mat_properties_[mb];

    for (int r = 0; r < (spec->regions).size(); r++) {
      std::string region = (spec->regions)[r];
      auto block = mesh_->getSetEntities(region, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

      AmanziMesh::Entity_ID_List::iterator c;
      if (phase == TRANSPORT_PHASE_LIQUID) {
        for (const auto& c: block) {
          D_[c] +=
            md * spec->tau[phase] * porosity[0][c] * saturation[0][c] * mol_density[0][c];
        }
      } else if (phase == TRANSPORT_PHASE_GAS) {
        for (const auto& c: block) {
          D_[c] += md * spec->tau[phase] * porosity[0][c] * (1.0 - saturation[0][c]) *
                    mol_density[0][c];
        }
      }
    }
  }
}


/* ******************************************************************
* Check all phases for the given name.
****************************************************************** */
int
Transport_ATS::FindDiffusionValue(const std::string& tcc_name, double* md, int* phase)
{
  for (int i = 0; i < TRANSPORT_NUMBER_PHASES; i++) {
    if (diffusion_phase_[i] == Teuchos::null) continue;
    int ok = diffusion_phase_[i]->FindDiffusionValue(tcc_name, md);
    if (ok == 0) {
      *phase = i;
      return 0;
    }
  }

  *md = 0.0;
  *phase = -1;
  return -1;
}


/* ******************************************************************
*  Find direction of axi-symmetry.
****************************************************************** */
void
Transport_ATS::CalculateAxiSymmetryDirection()
{
  axi_symmetry_.resize(ncells_owned, -1);
  if (S_->HasRecord(permeability_key_, tag_next_)) {
    const Epetra_MultiVector& perm =
      *S_->Get<CompositeVector>(permeability_key_, tag_next_).ViewComponent("cell");

    if (perm.NumVectors() == 3) {
      for (int c = 0; c < ncells_owned; ++c) {
        int k;
        if (perm[0][c] != perm[1][c] && perm[1][c] == perm[2][c]) {
          k = 0;
        } else if (perm[1][c] != perm[2][c] && perm[2][c] == perm[0][c]) {
          k = 1;
        } else if (perm[2][c] != perm[0][c] && perm[0][c] == perm[1][c]) {
          k = 2;
        }
        axi_symmetry_[c] = k;
      }
    }
  }
}

} // namespace Transport
} // namespace Amanzi
