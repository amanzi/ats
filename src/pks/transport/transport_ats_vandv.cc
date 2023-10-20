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

/* ****************************************************************
* Construct default state for unit tests.
**************************************************************** */
void
Transport_ATS::CreateDefaultState(Teuchos::RCP<const AmanziMesh::Mesh>& mesh, int ncomponents)
{
  std::string name("state");
  S_->Require<double>("fluid_density", Tags::DEFAULT, name);

  if (!S_->HasRecord(saturation_key_, Tags::NEXT)) {
    S_->Require<CompositeVector, CompositeVectorSpace>(saturation_key_, Tags::NEXT, name)
      .SetMesh(mesh)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  }

  if (!S_->HasRecord(saturation_key_, Tags::CURRENT)) {
    S_->Require<CompositeVector, CompositeVectorSpace>(saturation_key_, Tags::CURRENT, name)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  }

  if (!S_->HasRecord(flux_key_, Tags::NEXT)) {
    S_->Require<CompositeVector, CompositeVectorSpace>(flux_key_, Tags::NEXT, name)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->SetComponent("face", AmanziMesh::Entity_kind::FACE, 1);
  }

  if (!S_->HasRecord(tcc_key_, Tags::NEXT)) {
    std::vector<std::string> subfield_names;
    for (int i = 0; i != ncomponents; ++i) { subfield_names.push_back(component_names_[i]); }
    S_->Require<CompositeVector, CompositeVectorSpace>(tcc_key_, Tags::NEXT, name)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, ncomponents);
    S_->GetRecordSetW(tcc_key_).set_subfieldnames(subfield_names);
  }

  // initialize fields
  S_->Setup();

  // set popular default values
  S_->GetW<double>("fluid_density", Tags::DEFAULT, name) = 1000.0;
  S_->GetRecordW("fluid_density", Tags::DEFAULT, name).set_initialized();

  S_->GetW<CompositeVector>(saturation_key_, Tags::NEXT, name).PutScalar(1.0);
  S_->GetRecordW(saturation_key_, Tags::NEXT, name).set_initialized();

  S_->GetW<CompositeVector>(saturation_key_, Tags::CURRENT, name).PutScalar(1.0);
  S_->GetRecordW(saturation_key_, Tags::CURRENT, name).set_initialized();

  S_->GetW<CompositeVector>(tcc_key_, Tags::NEXT, name).PutScalar(0.0);
  S_->GetRecordW(tcc_key_, Tags::NEXT, name).set_initialized();

  S_->GetW<CompositeVector>(flux_key_, Tags::NEXT, name).PutScalar(0.0);
  S_->GetRecordW(flux_key_, Tags::NEXT, name).set_initialized();

  S_->InitializeFields();
}


/* *******************************************************************
* Routine verifies that the velocity field is divergence free
******************************************************************* */
void
Transport_ATS::Policy(const Tag& tag)
{
  if (mesh_->getComm()->NumProc() > 1) {
    if (!S_->Get<CompositeVector>(tcc_key_, tag).Ghosted()) {
      Errors::Message msg;
      msg << "Field \"total component concentration\" has no ghost values."
          << " Transport PK is giving up.\n";
      Exceptions::amanzi_throw(msg);
    }
  }
}


/* *******************************************************************
* Calculates extrema of specified solutes and print them.
******************************************************************* */
void
Transport_ATS::VV_PrintSoluteExtrema(const Epetra_MultiVector& tcc_next, double dT_MPC)
{
  int num_components = tcc_next.NumVectors();
  double tccmin_vec[num_components];
  double tccmax_vec[num_components];

  tcc_next.MinValue(tccmin_vec);
  tcc_next.MaxValue(tccmax_vec);

  for (int n = 0; n < runtime_solutes_.size(); n++) {
    int i = FindComponentNumber(runtime_solutes_[n]);
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

          auto cells = mesh_->getFaceCells(f, Amanzi::AmanziMesh::Parallel_kind::ALL);
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

    // old capability
    //mass_solutes_exact_[i] += VV_SoluteVolumeChangePerSecond(i) * dT_MPC;
    double mass_solute(0.0);
    for (int c = 0; c < ncells_owned; c++) {
      double vol = mesh_->getCellVolume(c);
      mass_solute += (*ws_)[0][c] * (*phi_)[0][c] * tcc_next[i][c] * vol * (*mol_dens_)[0][c];
    }
    // mass_solute /= units_.concentration_factor();
    // mass_solutes_stepstart_[i] /= units_.concentration_factor();
    // mass_solutes_bc_[i] /= units_.concentration_factor();

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
Transport_ATS::VV_CheckInfluxBC() const
{
  int number_components = tcc->ViewComponent("cell")->NumVectors();
  std::vector<int> influx_face(nfaces_wghost);

  for (int i = 0; i < number_components; i++) {
    influx_face.assign(nfaces_wghost, 0);

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
Transport_ATS::VV_CheckGEDproperty(Epetra_MultiVector& tracer) const
{
  int i, num_components = tracer.NumVectors();
  double tr_min[num_components];
  double tr_max[num_components];

  tracer.MinValue(tr_min);
  tracer.MaxValue(tr_max);

  for (i = 0; i < num_components; i++) {
    if (tr_min[i] < 0) {
      std::cout << "Transport_ATS: concentration violates GED property" << std::endl;
      std::cout << "    Make an Amanzi ticket or turn off internal transport tests" << std::endl;
      std::cout << "    MyPID = " << MyPID << std::endl;
      std::cout << "    component = " << i << std::endl;
      std::cout << "    time = " << t_physics_ << std::endl;
      std::cout << "    min/max values = " << tr_min[i] << " " << tr_max[i] << std::endl;

      Errors::Message msg;
      msg << "Concentration violates GED property."
          << "\n";
      Exceptions::amanzi_throw(msg);
    }
  }
}


/* *******************************************************************
 * Check that the tracer is between 0 and 1.
 ****************************************************************** */
void
Transport_ATS::VV_CheckTracerBounds(Epetra_MultiVector& tracer,
                                    int component,
                                    double lower_bound,
                                    double upper_bound,
                                    double tol) const
{
  Epetra_MultiVector& tcc_prev = *tcc->ViewComponent("cell");

  for (int c = 0; c < ncells_owned; c++) {
    double value = tracer[component][c];
    if (value < lower_bound - tol || value > upper_bound + tol) {
      std::cout << "Transport_ATS: tracer violates bounds" << std::endl;
      std::cout << "    Make an Amanzi ticket or turn off internal transport tests" << std::endl;
      std::cout << "    MyPID = " << MyPID << std::endl;
      std::cout << "    component = " << component << std::endl;
      std::cout << "    simulation time = " << t_physics_ << std::endl;
      std::cout << "      cell = " << c << std::endl;
      std::cout << "      center = " << mesh_->getCellCentroid(c) << std::endl;
      std::cout << "      value (old) = " << tcc_prev[component][c] << std::endl;
      std::cout << "      value (new) = " << value << std::endl;

      Errors::Message msg;
      msg << "Tracer violates bounds."
          << "\n";
      Exceptions::amanzi_throw(msg);
    }
  }
}


/* ******************************************************************
* Calculate change of tracer volume per second due to boundary flux.
* This is the simplified version (lipnikov@lanl.gov).
****************************************************************** */
double
Transport_ATS::VV_SoluteVolumeChangePerSecond(int idx_tracer)
{
  double volume = 0.0;

  for (int m = 0; m < bcs_.size(); m++) {
    std::vector<int>& tcc_index = bcs_[m]->tcc_index();
    int ncomp = tcc_index.size();

    for (int i = 0; i < ncomp; i++) {
      if (tcc_index[i] == idx_tracer) {
        for (auto it = bcs_[m]->begin(); it != bcs_[m]->end(); ++it) {
          int f = it->first;
          std::vector<double>& values = it->second;

          int c2 = (*downwind_cell_)[f];

          if (f < nfaces_owned && c2 >= 0) {
            double u = fabs((*flux_)[0][f]);
            volume += u * values[i];
          }
        }
      }
    }
  }
  return volume;
}

/* *******************************************************************
* Error estimate uses analytic function and solution.
* ***************************************************************** */
void
Transport_ATS::CalculateLpErrors(AnalyticFunction f,
                                 double t,
                                 Epetra_Vector* sol,
                                 double* L1,
                                 double* L2)
{
  *L1 = *L2 = 0.0;
  for (int c = 0; c < sol->MyLength(); c++) {
    const AmanziGeometry::Point& xc = mesh_->getCellCentroid(c);
    double d = (*sol)[c] - f(xc, t);

    double volume = mesh_->getCellVolume(c);
    *L1 += fabs(d) * volume;
    *L2 += d * d * volume;
  }

  *L2 = sqrt(*L2);
}


double
Transport_ATS::ComputeSolute(const Epetra_MultiVector& tcc_c, int i)
{
  double mass_solute(0.0);
  for (int c = 0; c < ncells_owned; c++) {
    double vol = mesh_->getCellVolume(c);
    mass_solute += (*ws_next)[0][c] * (*phi_)[0][c] * tcc_c[i][c] * vol * (*mol_dens_next)[0][c] +
                   (*solid_qty_)[i][c];
  }
  //mass_solute /= units_.concentration_factor();

  double tmp1 = mass_solute;
  mesh_->getComm()->SumAll(&tmp1, &mass_solute, 1);

  return mass_solute;
}


double
Transport_ATS::ComputeSolute(const Epetra_MultiVector& tcc_c,
                             const Epetra_MultiVector& ws,
                             const Epetra_MultiVector& den,
                             int i)
{
  double mass_solute(0.0);
  for (int c = 0; c < ncells_owned; c++) {
    double vol = mesh_->getCellVolume(c);
    mass_solute += ws[0][c] * (*phi_)[0][c] * tcc_c[i][c] * vol * den[0][c] + (*solid_qty_)[i][c];
  }
  //mass_solute /= units_.concentration_factor();

  double tmp1 = mass_solute;
  mesh_->getComm()->SumAll(&tmp1, &mass_solute, 1);

  return mass_solute;
}

} // namespace Transport
} // namespace Amanzi
