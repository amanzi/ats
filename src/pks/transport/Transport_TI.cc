/*
  Transport PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <algorithm>

#include "ReconstructionCell.hh"
#include "OperatorDefs.hh"
#include "Transport_PK_ATS.hh"

namespace Amanzi {
namespace Transport {

/* ******************************************************************* 
 * Routine takes a parallel overlapping vector C and returns a parallel
 * overlapping vector F(C).
 ****************************************************************** */

void Transport_PK_ATS::FunctionalTimeDerivative(double t,
                                                const Epetra_Vector& component,
                                                Epetra_Vector& f_component)

{

  // distribute vector
  auto component_tmp = Teuchos::rcp(new Epetra_Vector(component));
  component_tmp->Import(component, *tcc->importer("cell"), Insert);
  
  Teuchos::ParameterList plist = tp_list_->sublist("reconstruction");
  lifting_->Init(plist);
  lifting_->ComputeGradient(component_tmp);

  // extract boundary conditions for the current component
  std::vector<int> bc_model(nfaces_wghost, Operators::OPERATOR_BC_NONE);
  std::vector<double> bc_value(nfaces_wghost);

  for (int m = 0; m < bcs_.size(); m++) {
    std::vector<int>& tcc_index = bcs_[m]->tcc_index();
    int ncomp = tcc_index.size();

    for (int i = 0; i < ncomp; i++) {
      if (current_component_ == tcc_index[i]) {
        for (auto it = bcs_[m]->begin(); it != bcs_[m]->end(); ++it) {
          int f = it->first;
          std::vector<double>& values = it->second;

          bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
          bc_value[f] = values[i];
        }
      }
    }
  }


  limiter_->Init(plist, flux_);
  limiter_->ApplyLimiter(component_tmp, 0, lifting_->gradient(), bc_model, bc_value); 
  limiter_->gradient()->ScatterMasterToGhosted("cell");
  

  // ADVECTIVE FLUXES
  // We assume that limiters made their job up to round-off errors.
  // Min-max condition will enforce robustness w.r.t. these errors.
  int f, c1, c2;
  double u, u1, u2, umin, umax, upwind_tcc, tcc_flux;

  f_component.PutScalar(0.0);
  for (int f = 0; f < nfaces_wghost; f++) {  // loop over master and slave faces
    c1 = (*upwind_cell_)[f];
    c2 = (*downwind_cell_)[f];
  

    if (c1 >= 0 && c2 >= 0) {
      u1 = component[c1];
      u2 = component[c2];
      umin = std::min(u1, u2);
      umax = std::max(u1, u2);
    } else if (c1 >= 0) {
      u1 = u2 = umin = umax = component[c1];
    } else if (c2 >= 0) {
      u1 = u2 = umin = umax = component[c2];
    }

    u = fabs((*flux_)[0][f]);
    const AmanziGeometry::Point& xf = mesh_->face_centroid(f);

    if (c1 >= 0 && c1 < ncells_owned && c2 >= 0 && c2 < ncells_owned) {
      upwind_tcc = limiter_->getValue(c1, xf);
      upwind_tcc = std::max(upwind_tcc, umin);
      upwind_tcc = std::min(upwind_tcc, umax);

      tcc_flux = u * upwind_tcc;
      f_component[c1] -= tcc_flux;
      f_component[c2] += tcc_flux;

    } else if (c1 >= 0 && c1 < ncells_owned && (c2 >= ncells_owned || c2 < 0)) {

      upwind_tcc = limiter_->getValue(c1, xf);
      upwind_tcc = std::max(upwind_tcc, umin);
      upwind_tcc = std::min(upwind_tcc, umax);

      tcc_flux = u * upwind_tcc;
      f_component[c1] -= tcc_flux;
    } else if (c1 >= 0 && c1 < ncells_owned && (c2 < 0)) {
      upwind_tcc = component[c1];
      upwind_tcc = std::max(upwind_tcc, umin);
      upwind_tcc = std::min(upwind_tcc, umax);

      tcc_flux = u * upwind_tcc;
      f_component[c1] -= tcc_flux;
      
    } else if (c1 >= ncells_owned && c2 >= 0 && c2 < ncells_owned) {
      upwind_tcc = limiter_->getValue(c1, xf);
      upwind_tcc = std::max(upwind_tcc, umin);
      upwind_tcc = std::min(upwind_tcc, umax);

      tcc_flux = u * upwind_tcc;
      f_component[c2] += tcc_flux;
    }
  }

                                                

  // process external sources
  if (srcs_.size() != 0) {
    ComputeAddSourceTerms(t, 1., f_component, current_component_, current_component_);
  }

  // if (domain_ == "surface"){
  //   double mass = 0., mass_mov_loc = 0;
  //   for (int c = 0; c < ncells_owned; c++)  mass_mov_loc += f_component[c];
  //   mesh_->get_comm()->SumAll(&mass_mov_loc, &mass, 1);
  //   *vo_->os()<<"mass mov "<< mass<<"\n";
  // }

  for (int c = 0; c < ncells_owned; c++) {  // calculate conservative quantatity   
    double vol_phi_ws_den = mesh_->cell_volume(c) * (*phi_)[0][c] * (*ws_start)[0][c] * (*mol_dens_start)[0][c];
    if ((*ws_start)[0][c] < 1e-12)
      vol_phi_ws_den = mesh_->cell_volume(c) * (*phi_)[0][c] * (*ws_end)[0][c] * (*mol_dens_end)[0][c];

    if (vol_phi_ws_den > water_tolerance_){
      f_component[c] /= vol_phi_ws_den;
    }
  }

  // BOUNDARY CONDITIONS for ADVECTION
  for (int m = 0; m < bcs_.size(); m++) {
    std::vector<int>& tcc_index = bcs_[m]->tcc_index();
    int ncomp = tcc_index.size();

    for (int i = 0; i < ncomp; i++) {
      if (current_component_ == tcc_index[i]) {
        for (auto it = bcs_[m]->begin(); it != bcs_[m]->end(); ++it) {
          int f = it->first;
          std::vector<double>& values = it->second;
          c2 = (*downwind_cell_)[f];

          if (c2 >= 0 && f < nfaces_owned) {
            u = fabs((*flux_)[0][f]);
            double vol_phi_ws_den = mesh_->cell_volume(c2) * (*phi_)[0][c2] * (*ws_start)[0][c2] * (*mol_dens_start)[0][c2];
            if ((*ws_start)[0][c2] < 1e-12)
              vol_phi_ws_den = mesh_->cell_volume(c2) * (*phi_)[0][c2] * (*ws_end)[0][c2] * (*mol_dens_end)[0][c2];

            tcc_flux = u * values[i];

            //if (abs(tcc_flux) > 1e-9)
            // if ((mesh_->face_centroid(f)[0]>9.5)&&(mesh_->face_centroid(f)[0]<10.5))              
            //   std::cout<<mesh_->get_comm()->MyPID()<<" c2 "<<c2<<" mass "<<tcc_flux<<" "<<mesh_->face_centroid(f)<<"\n";
            
            if (vol_phi_ws_den > water_tolerance_ ){
              f_component[c2] += tcc_flux / vol_phi_ws_den;            
            }
            
          }
        }
      }
    }
  }
}


  
}  // namespace Transport
}  // namespace Amanzi


