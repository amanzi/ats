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

#include <algorithm>
#include <vector>

#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "errors.hh"

#include "transport_ats.hh"

namespace Amanzi {
namespace Transport {

/* *******************************************************************
* Calculates extrema of specified solutes and print them.
******************************************************************* */
void
Transport_ATS::PrintSoluteExtrema(const Epetra_MultiVector& tcc_next, double dT_MPC)
{
  int num_components = tcc_next.NumVectors();
  double tccmin_vec[num_components];
  double tccmax_vec[num_components];

  tcc_next.MinValue(tccmin_vec);
  tcc_next.MaxValue(tccmax_vec);

  for (int n = 0; n < runtime_solutes_.size(); n++) {
    int i = FindComponentNumber_(runtime_solutes_[n]);
    double tccmin, tccmax;
    tcc_next.Comm().MinAll(&(tccmin_vec[i]), &tccmin, 1); // find the global extrema
    tcc_next.Comm().MaxAll(&(tccmax_vec[i]), &tccmax, 1);

    int nregions = runtime_regions_.size();
    double solute_flux(0.0);
    bool flag(false);

    for (int k = 0; k < nregions; k++) {
      if (mesh_->isValidSetName(runtime_regions_[k], AmanziMesh::Entity_kind::FACE)) {
        flag = true;
        auto block = mesh_->getSetEntities(
          runtime_regions_[k], AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
        int nblock = block.size();

        for (int m = 0; m < nblock; m++) {
          int f = block[m];

          auto cells = mesh_->getFaceCells(f);
          int dir, c = cells[0];

          const AmanziGeometry::Point& normal = mesh_->getFaceNormal(f, c, &dir);
          double u = (*flux_)[0][f] * dir;
          if (u > 0) solute_flux += u * tcc_next[i][c];
        }
      }
    }

    //solute_flux *= units_.concentration_factor();

    double tmp = solute_flux;
    mesh_->getComm()->SumAll(&tmp, &solute_flux, 1);

    double ws_min, ws_max;
    ws_->MinValue(&ws_min);
    ws_->MaxValue(&ws_max);

    *vo_->os() << runtime_solutes_[n] << ": min=" << tccmin << " max=" << tccmax << " ws: "
               << "min=" << ws_min << " max=" << ws_max << "\n";
    if (flag) *vo_->os() << ", flux=" << solute_flux << " mol/s";

    double mass_solute(0.0);
    for (int c = 0; c < tcc_next.MyLength(); c++) {
      double vol = mesh_->getCellVolume(c);
      mass_solute += (*ws_)[0][c] * (*phi_)[0][c] * tcc_next[i][c] * vol * (*mol_dens_)[0][c];
    }

    double tmp1 = mass_solute, tmp2 = mass_solutes_exact_[i], mass_exact;
    double tmp_start = mass_solutes_stepstart_[i];
    double tmp_bc = mass_solutes_bc_[i];
    mesh_->getComm()->SumAll(&tmp1, &mass_solute, 1);
    mesh_->getComm()->SumAll(&tmp2, &mass_exact, 1);
    mesh_->getComm()->SumAll(&tmp_start, &(mass_solutes_stepstart_[i]), 1);
    mesh_->getComm()->SumAll(&tmp_bc, &(mass_solutes_bc_[i]), 1);
  }
}


/********************************************************************
* Check completeness of influx boundary conditions.
****************************************************************** */
void
Transport_ATS::CheckInfluxBC_() const
{
  int nfaces_all =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::ALL);
  int number_components = tcc->ViewComponent("cell")->NumVectors();
  std::vector<int> influx_face(nfaces_all);

  for (int i = 0; i < number_components; i++) {
    influx_face.assign(nfaces_all, 0);

    for (int m = 0; m < bcs_.size(); m++) {
      std::vector<int>& tcc_index = bcs_[m]->tcc_index();
      int ncomp = tcc_index.size();

      for (int k = 0; k < ncomp; k++) {
        if (i == tcc_index[k]) {
          for (auto it = bcs_[m]->begin(); it != bcs_[m]->end(); ++it) {
            int f = it->first;
            influx_face[f] = 1;
          }
        }
      }
    }

    for (int m = 0; m < bcs_.size(); m++) {
      std::vector<int>& tcc_index = bcs_[m]->tcc_index();
      int ncomp = tcc_index.size();

      for (int k = 0; k < ncomp; k++) {
        if (i == tcc_index[k]) {
          for (auto it = bcs_[m]->begin(); it != bcs_[m]->end(); ++it) {
            int f = it->first;
            if ((*flux_)[0][f] < 0 && influx_face[f] == 0) {
              Errors::Message msg;
              msg << "No influx boundary condition has been found for component " << i << ".\n";
              Exceptions::amanzi_throw(msg);
            }
          }
        }
      }
    }
  }
}


/* *******************************************************************
 * Check that global extrema diminished
 ****************************************************************** */
void
CheckGEDProperty(const Epetra_MultiVector& tracer, double t_physics)
{
  int i, num_components = tracer.NumVectors();
  double tr_min[num_components];

  tracer.MinValue(tr_min);
  double tr_all_min = *std::min_element(&tr_min[0], &tr_min[num_components]);

  if (tr_all_min < 0) {
    double tr_max[num_components];
    tracer.MaxValue(tr_max);
    int MyPID = tracer.Comm().MyPID();
    for (i = 0; i < num_components; i++) {
      if (tr_min[i] < 0) {
        std::cout << "Transport_ATS: concentration violates GED property" << std::endl;
        std::cout << "    Make an Amanzi ticket or turn off internal transport tests" << std::endl;
        std::cout << "    MyPID = " << MyPID << std::endl;
        std::cout << "    component = " << i << std::endl;
        std::cout << "    time = " << t_physics << std::endl;
        std::cout << "    min/max values = " << tr_min[i] << " " << tr_max[i] << std::endl;

        Errors::Message msg;
        msg << "Concentration violates GED property."
            << "\n";
        Exceptions::amanzi_throw(msg);
      }
    }
  }
}


/* *******************************************************************
 * Check that the tracer is between 0 and 1.
 ****************************************************************** */
void
CheckTracerBounds(const Epetra_MultiVector& tcc,
                  const Epetra_MultiVector& tcc_prev,
                  const AmanziMesh::Mesh& mesh,
                  double t_physics,
                  int component,
                  double lower_bound,
                  double upper_bound,
                  double tol)
{
  for (int c = 0; c < tcc_prev.MyLength(); c++) {
    double value = tcc[component][c];
    if (value < lower_bound - tol || value > upper_bound + tol) {
      int MyPID = tcc.Comm().MyPID();
      std::cout << "Transport_ATS: TCC violates bounds" << std::endl;
      std::cout << "    Make an Amanzi ticket or turn off internal transport tests" << std::endl;
      std::cout << "    MyPID = " << MyPID << std::endl;
      std::cout << "    component = " << component << std::endl;
      std::cout << "    simulation time = " << t_physics << std::endl;
      std::cout << "      cell = " << c << std::endl;
      std::cout << "      center = " << mesh.getCellCentroid(c) << std::endl;
      std::cout << "      value (old) = " << tcc_prev[component][c] << std::endl;
      std::cout << "      value (new) = " << value << std::endl;

      Errors::Message msg;
      msg << "TCC violates bounds."
          << "\n";
      Exceptions::amanzi_throw(msg);
    }
  }
}

} // namespace Transport
} // namespace Amanzi
