/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  An elevation evaluator getting values from the columnar meshes.

  Authors: Ahmad Jan (jana@ornl.gov)
*/
#include <boost/algorithm/string/predicate.hpp>

#include "Mesh.hh"
#include "Point.hh"
#include "elevation_evaluator_column.hh"

namespace Amanzi {
namespace Flow {

ElevationEvaluatorColumn::ElevationEvaluatorColumn(Teuchos::ParameterList& plist) :
  ElevationEvaluator(plist) 
{
  surface_domain_ =  plist.get<std::string>("surface domain name", "surface_star");
  dset_name_ = plist.get<std::string>("domain set name", "column");
};

ElevationEvaluatorColumn::ElevationEvaluatorColumn(const ElevationEvaluatorColumn& other) :
  ElevationEvaluator(other),
  base_por_key_(other.base_por_key_),
  surface_domain_(other.surface_domain_),
  dset_name_(other.dset_name_)
{};

Teuchos::RCP<FieldEvaluator>
ElevationEvaluatorColumn::Clone() const {
  return Teuchos::rcp(new ElevationEvaluatorColumn(*this));
}

void ElevationEvaluatorColumn::EvaluateElevationAndSlope_(const Teuchos::Ptr<State>& S,
        const std::vector<Teuchos::Ptr<CompositeVector> >& results) {

  Teuchos::Ptr<CompositeVector> elev = results[0];
  Teuchos::Ptr<CompositeVector> slope = results[1];
  Epetra_MultiVector& elev_c = *elev->ViewComponent("cell", false);
  Epetra_MultiVector& slope_c = *slope->ViewComponent("cell", false);
 
  // Get the elevation and slope values from the domain mesh.
  Key domain_sf = Keys::getDomain(my_keys_[0]);
  
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  // Set the elevation on cells by getting the corresponding face and its
  // centroid.
  int ncells = elev_c.MyLength();
  std::vector<AmanziGeometry::Point> my_centroid;

  auto domain_set = S->GetDomainSet(dset_name_);
  const Epetra_Map& cell_map = S->GetMesh(surface_domain_)->cell_map(false);
  
  for (const auto& domain : *domain_set) {
    int id = Keys::getDomainSetIndex<int>(domain);
    int c = cell_map.LID(id);
    std::vector<AmanziGeometry::Point> coord;

    S->GetMesh(domain)->face_get_coordinates(0, &coord); // 0 is the id of top face of the column mesh
    
    elev_c[0][c] = coord[0][2];
  }
  
  //Now get slope
  if (domain_sf == surface_domain_){
    Teuchos::RCP<CompositeVector> elev_ngb = S->GetW<CompositeVector>(my_keys_[0], S->GetField(my_keys_[0]).owner() );
    elev_ngb->ScatterMasterToGhosted("cell");
    const Epetra_MultiVector& elev_ngb_c = *elev_ngb->ViewComponent("cell",true);
    

    //get all cell centroids
    for (int c=0; c!=ncells; ++c) {
      std::stringstream my_name;
      int id = S->GetMesh(surface_domain_)->cell_map(true).GID(c);
      AmanziGeometry::Point P1 = S->GetMesh(surface_domain_)->cell_centroid(c);
      P1.set(P1[0], P1[1], elev_ngb_c[0][c]);
      my_centroid.push_back(P1);
    }
    
    //get neighboring cell ids
    for (int c=0; c!= ncells; c++){
      AmanziMesh::Entity_ID_List nadj_cellids;
      S->GetMesh(surface_domain_)->cell_get_face_adj_cells(c, AmanziMesh::Parallel_type::ALL, &nadj_cellids);
      int nface_pcell = S->GetMesh(surface_domain_)->cell_get_num_faces(c);

      int ngb_cells = nadj_cellids.size();
      std::vector<AmanziGeometry::Point> ngb_centroids(ngb_cells);
          
      //get the neighboring cell's centroids
      for(unsigned i=0; i<ngb_cells; i++){
        AmanziGeometry::Point P2 = S->GetMesh(surface_domain_)->cell_centroid(nadj_cellids[i]);
        ngb_centroids[i].set(P2[0], P2[1], elev_ngb_c[0][nadj_cellids[i]]);
      }
      
      int id = S->GetMesh(surface_domain_)->cell_map(false).GID(c);
      Key my_name = Keys::getDomainInSet(dset_name_, id);
      
      std::vector<AmanziGeometry::Point> Normal;
      AmanziGeometry::Point N, PQ, PR, Nor_avg(3);
      
      if (ngb_cells >1){
	for (int i=0; i <ngb_cells-1; i++){
	  PQ = my_centroid[c] - ngb_centroids[i];
	  PR = my_centroid[c] - ngb_centroids[i+1];
	  N = PQ^PR;
          if (N[2] < 0)
            N *= -1.; // all normals upward
	  Normal.push_back(N);
	}

        AmanziGeometry::Point fnor = S->GetMesh(my_name)->face_normal(0); //0 is the id of top face
	Nor_avg = (nface_pcell - Normal.size()) * fnor; 
	for (int i=0; i <Normal.size(); i++)
	  Nor_avg += Normal[i];
	
        Nor_avg /= nface_pcell;
	slope_c[0][c] = (std::sqrt(std::pow(Nor_avg[0],2) + std::pow(Nor_avg[1],2)))/ std::abs(Nor_avg[2]);
        
      }
      else if (ngb_cells == 1){
	PQ = my_centroid[c] - ngb_centroids[0];
	slope_c[0][c] = std::abs(PQ[2]) / (std::sqrt(std::pow(PQ[0],2) + std::pow(PQ[1],2))); 
      }
      else if (ngb_cells == 0){
	slope_c[0][c] = 0.0;
      }
    }
    
  }
  else{
    slope_c[0][0] = 0.0; // if domain is surface_column:*, slope is zero.
  }

  
  if (elev->HasComponent("face")) {
    Epetra_MultiVector& elev_f = *elev->ViewComponent("face", false);

    Teuchos::RCP<CompositeVector> elev_ngb = S->GetW<CompositeVector>(my_keys_[0], S->GetField(my_keys_[0]).owner() );
    elev_ngb->ScatterMasterToGhosted("cell");
    const Epetra_MultiVector& elev_ngb_c = *elev_ngb->ViewComponent("cell",true);
    int nfaces = elev_f.MyLength();
    
    for (int f=0; f!=nfaces; ++f) {
      AmanziMesh::Entity_ID_List nadj_cellids;
      S->GetMesh(surface_domain_)->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &nadj_cellids);
      double ef = 0;
      for (int i=0; i<nadj_cellids.size(); i++){
        ef += elev_ngb_c[0][nadj_cellids[i]];
      }
      elev_f[0][f] = ef/nadj_cellids.size();
    }
  }
  
}


void ElevationEvaluatorColumn::EnsureCompatibility(const Teuchos::Ptr<State>& S){
  
  Key domain_sf = Keys::getDomain(my_keys_[0]);

  int ncells = S->GetMesh(surface_domain_)->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
 
  if (domain_sf == surface_domain_) {
    for (int c =0; c < ncells; c++){
      int id = S->GetMesh(surface_domain_)->cell_map(false).GID(c);
      Key name = Keys::getDomainInSet(dset_name_, id);
      base_por_key_ = Keys::readKey(plist_, name, "base porosity", "base_porosity");
      dependencies_.insert(base_por_key_);
    }
  } else {
    Errors::Message msg("ElevationEvaluatorColumn: this evaluator should be used for columnar meshes only.");
    Exceptions::amanzi_throw(msg);
  }

  ElevationEvaluator::EnsureCompatibility(S.ptr());
}
} //namespace
} //namespace
