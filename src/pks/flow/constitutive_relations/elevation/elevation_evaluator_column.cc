/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ahmad Jan (jana@ornl.gov)
*/

/*
  An elevation evaluator getting values from the columnar meshes.

*/

#include "Mesh.hh"
#include "Point.hh"
#include "elevation_evaluator_column.hh"

namespace Amanzi {
namespace Flow {

ColumnElevationEvaluator::ColumnElevationEvaluator(Teuchos::ParameterList& plist)
  : ElevationEvaluator(plist)
{
  dset_name_ = plist.get<std::string>("domain set name", "column");
  surface_domain_ = Keys::getDomain(my_keys_.front().first);
  base_poro_suffix_ = Keys::readSuffix(plist_, "base porosity", "base_porosity");
};


Teuchos::RCP<Evaluator>
ColumnElevationEvaluator::Clone() const
{
  return Teuchos::rcp(new ColumnElevationEvaluator(*this));
}


void
ColumnElevationEvaluator::EnsureEvaluators(State& S)
{
  Tag tag = Keys::readTag(plist_, my_keys_.front().second);
  auto dset = S.GetDomainSet(dset_name_);
  for (const auto& domain : *dset) {
    dependencies_.insert(KeyTag{ Keys::getKey(domain, base_poro_suffix_), tag });
  }
  ElevationEvaluator::EnsureEvaluators(S);
}


void
ColumnElevationEvaluator::EvaluateElevationAndSlope_(const State& S,
                                                     const std::vector<CompositeVector*>& results)
{
  CompositeVector* elev = results[0];
  CompositeVector* slope = results[1];
  CompositeVector* aspect = results[2];

  aspect->PutScalar(0.); // no aspect to a column mesh?
  Epetra_MultiVector& elev_c = *elev->ViewComponent("cell", false);
  Epetra_MultiVector& slope_c = *slope->ViewComponent("cell", false);

  // Set the elevation on cells by getting the corresponding face and its
  // centroid.
  int ncells = elev_c.MyLength();
  std::vector<AmanziGeometry::Point> my_centroid;

  auto domain_set = S.GetDomainSet(dset_name_);
  const Epetra_Map& cell_map = S.GetMesh(surface_domain_)->getMap(AmanziMesh::Entity_kind::CELL,false);

  AMANZI_ASSERT(domain_set->size() == cell_map.NumMyElements());
  AMANZI_ASSERT(domain_set->size() == elev_c.MyLength());
  for (const auto& domain : *domain_set) {
    int gid = Keys::getDomainSetIndex<int>(domain);
    int c = cell_map.LID(gid);
    auto coord = S.GetMesh(domain)->getFaceCentroid(0); // 0 is the id of top face of the column mesh
    elev_c[0][c] = coord[2];
  }

  // Now get slope
  elev->ScatterMasterToGhosted("cell");
  const Epetra_MultiVector& elev_ngb_c = *elev->ViewComponent("cell", true);

  // get all cell centroids
  for (int c = 0; c != ncells; ++c) {
    int id = S.GetMesh(surface_domain_)->getMap(AmanziMesh::Entity_kind::CELL,true).GID(c);
    AmanziGeometry::Point P1 = S.GetMesh(surface_domain_)->getCellCentroid(c);
    P1.set(P1[0], P1[1], elev_ngb_c[0][c]);
    my_centroid.push_back(P1);
  }

  // get neighboring cell ids
  for (int c = 0; c != ncells; c++) {
    auto nadj_cellids = AmanziMesh::MeshAlgorithms::getCellFaceAdjacentCells(*S.GetMesh(surface_domain_), c, AmanziMesh::Parallel_kind::ALL);
    int nface_pcell = S.GetMesh(surface_domain_)->getCellNumFaces(c);

    int ngb_cells = nadj_cellids.size();
    std::vector<AmanziGeometry::Point> ngb_centroids(ngb_cells);

    // get the neighboring cell's centroids
    for (unsigned i = 0; i < ngb_cells; i++) {
      AmanziGeometry::Point P2 = S.GetMesh(surface_domain_)->getCellCentroid(nadj_cellids[i]);
      ngb_centroids[i].set(P2[0], P2[1], elev_ngb_c[0][nadj_cellids[i]]);
    }

    int id = S.GetMesh(surface_domain_)->getMap(AmanziMesh::Entity_kind::CELL,false).GID(c);
    Key my_name = Keys::getDomainInSet(dset_name_, id);

    std::vector<AmanziGeometry::Point> Normal;
    AmanziGeometry::Point N, PQ, PR, Nor_avg(3);

    if (ngb_cells > 1) {
      for (int i = 0; i < ngb_cells - 1; i++) {
        PQ = my_centroid[c] - ngb_centroids[i];
        PR = my_centroid[c] - ngb_centroids[i + 1];
        N = PQ ^ PR;
        if (N[2] < 0) N *= -1.; // all normals upward
        Normal.push_back(N);
      }

      AmanziGeometry::Point fnor = S.GetMesh(my_name)->getFaceNormal(0); //0 is the id of top face
      Nor_avg = (nface_pcell - Normal.size()) * fnor;
      for (int i = 0; i < Normal.size(); i++) Nor_avg += Normal[i];

      Nor_avg /= nface_pcell;
      slope_c[0][c] =
        (std::sqrt(std::pow(Nor_avg[0], 2) + std::pow(Nor_avg[1], 2))) / std::abs(Nor_avg[2]);

    } else if (ngb_cells == 1) {
      PQ = my_centroid[c] - ngb_centroids[0];
      slope_c[0][c] = std::abs(PQ[2]) / (std::sqrt(std::pow(PQ[0], 2) + std::pow(PQ[1], 2)));
    } else if (ngb_cells == 0) {
      slope_c[0][c] = 0.0;
    }
  }

  if (elev->HasComponent("face")) {
    Epetra_MultiVector& elev_f = *elev->ViewComponent("face", false);
    const Epetra_MultiVector& elev_ngb_c = *elev->ViewComponent("cell", true);
    int nfaces = elev_f.MyLength();

    for (int f = 0; f != nfaces; ++f) {
      auto nadj_cellids = S.GetMesh(surface_domain_)->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
      double ef = 0;
      for (int i = 0; i < nadj_cellids.size(); i++) { ef += elev_ngb_c[0][nadj_cellids[i]]; }
      elev_f[0][f] = ef / nadj_cellids.size();
    }
  }
}


} // namespace Flow
} // namespace Amanzi
