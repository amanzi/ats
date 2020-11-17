/*
  Transport PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_ATS_TRANSPORT_PK_HH_
#define AMANZI_ATS_TRANSPORT_PK_HH_

// TPLs
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Import.h"
#include "Teuchos_RCP.hpp"

// Amanzi
#include "CompositeVector.hh"
#include "DiffusionPhase.hh"
#include "Explicit_TI_FnBase.hh"
#include "MaterialProperties.hh"
#include "PK.hh"
#include "PK_Factory.hh"
#include "ReconstructionCell.hh"
#include "State.hh"
#include "Tensor.hh"
#include "Units.hh"
#include "VerboseObject.hh"
#include "Debugger.hh"
#include "PK_PhysicalExplicit.hh"
#include "DenseVector.hh"

#include <string>

#ifdef ALQUIMIA_ENABLED
#include "Alquimia_PK.hh"
#include "ChemistryEngine.hh"
#endif

// Transport
#include "LimiterCell.hh"
#include "MDMPartition.hh"
#include "MultiscaleTransportPorosityPartition.hh"
#include "TransportDomainFunction.hh"
#include "TransportDefs.hh"


/* ******************************************************************
The transport PK receives a reduced (optional) copy of a physical 
state at time n and returns a different state at time n+1. 

Unmodified physical quantaties in the returned state are the smart 
pointers to the original variables.
****************************************************************** */

namespace Amanzi {
namespace Transport {

typedef double AnalyticFunction(const AmanziGeometry::Point&, const double);

class Transport_ATS : public PK_PhysicalExplicit<Epetra_Vector> {

 public:
  Transport_ATS(Teuchos::ParameterList& pk_tree,
               const Teuchos::RCP<Teuchos::ParameterList>& glist,
               const Teuchos::RCP<State>& S,
               const Teuchos::RCP<TreeVector>& soln);

    Transport_ATS(const Teuchos::RCP<Teuchos::ParameterList>& glist,
                 Teuchos::RCP<State> S,
                 const std::string& pk_list_name,
                 std::vector<std::string>& component_names);

    ~Transport_ATS() = default;

  // members required by PK interface
  virtual void Setup(const Teuchos::Ptr<State>& S);
  virtual void Initialize(const Teuchos::Ptr<State>& S);

  virtual double get_dt();
  virtual void set_dt(double dt) {};

  virtual bool AdvanceStep(double t_old, double t_new, bool reinit=false); 
  virtual void CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S);
  virtual void CalculateDiagnostics(const Teuchos::RCP<State>& S) {};

  virtual void set_states(const Teuchos::RCP<State>& S,
                          const Teuchos::RCP<State>& S_inter,
                          const Teuchos::RCP<State>& S_next);

  virtual std::string name() { return "transport_ats"; }
  Key domain_name() {return domain_name_;}

  // main transport members
  // -- calculation of a stable time step needs saturations and darcy flux
  double StableTimeStep();
  void Sinks2TotalOutFlux(Epetra_MultiVector& tcc,
                          std::vector<double>& total_outflux, int n0, int n1);

  // coupling with chemistry
#ifdef ALQUIMIA_ENABLED
  void SetupAlquimia(Teuchos::RCP<AmanziChemistry::Alquimia_PK> chem_pk,
                     Teuchos::RCP<AmanziChemistry::ChemistryEngine> chem_engine);
#endif

  // -- access members  
  inline double cfl() { return cfl_; }
  Teuchos::RCP<const State> state() { return S_; }
  Teuchos::RCP<CompositeVector> total_component_concentration() { return tcc_tmp; }

  // -- control members
  void CreateDefaultState(Teuchos::RCP<const AmanziMesh::Mesh>& mesh, int ncomponents);
  void Policy(Teuchos::Ptr<State> S);

  void VV_CheckGEDproperty(Epetra_MultiVector& tracer) const; 
  void VV_CheckTracerBounds(Epetra_MultiVector& tracer, int component,
                            double lower_bound, double upper_bound, double tol = 0.0) const;
  void VV_CheckInfluxBC() const;
  void VV_PrintSoluteExtrema(const Epetra_MultiVector& tcc_next, double dT_MPC);
  double VV_SoluteVolumeChangePerSecond(int idx_solute);

  double ComputeSolute(const Epetra_MultiVector& tcc, int idx);

  double ComputeSolute(const Epetra_MultiVector& tcc,
                       const Epetra_MultiVector& ws,
                       const Epetra_MultiVector& den,
                       int idx);

  void CalculateLpErrors(AnalyticFunction f, double t, Epetra_Vector* sol, double* L1, double* L2);

  // -- sources and sinks for components from n0 to n1 including
  void ComputeAddSourceTerms(double tp, double dtp, 
                             Epetra_MultiVector& tcc, int n0, int n1);

  // void MixingSolutesWthSources(double told, double tnew);

  bool PopulateBoundaryData(std::vector<int>& bc_model,
                            std::vector<double>& bc_value, int component);

  // -- limiters 
  void LimiterBarthJespersen(const int component,
                             Teuchos::RCP<const Epetra_Vector> scalar_field, 
                             Teuchos::RCP<CompositeVector>& gradient, 
                             Teuchos::RCP<Epetra_Vector>& limiter);

  const std::vector<std::string>  component_names(){return component_names_;};
    int num_aqueous_component() {return num_aqueous;};
    int num_gaseous_component() {return num_gaseous;};    


 private:
  void InitializeFields_(const Teuchos::Ptr<State>& S);

  // advection members
  void AdvanceDonorUpwind(double dT);
  void AdvanceSecondOrderUpwindRKn(double dT);
  void AdvanceSecondOrderUpwindRK1(double dT);
  void AdvanceSecondOrderUpwindRK2(double dT);
  void Advance_Dispersion_Diffusion(double t_old, double t_new);

  // time integration members
  void FunctionalTimeDerivative(const double t, const Epetra_Vector& component, Epetra_Vector& f_component);
    //  void FunctionalTimeDerivative(const double t, const Epetra_Vector& component, TreeVector& f_component);

  void IdentifyUpwindCells();

  void InterpolateCellVector(
      const Epetra_MultiVector& v0, const Epetra_MultiVector& v1, 
      double dT_int, double dT, Epetra_MultiVector& v_int);

  const Teuchos::RCP<Epetra_IntVector>& upwind_cell() { return upwind_cell_; }
  const Teuchos::RCP<Epetra_IntVector>& downwind_cell() { return downwind_cell_; }  

  // physical models
  // -- dispersion and diffusion
  void CalculateDispersionTensor_(
      const Epetra_MultiVector& darcy_flux, const Epetra_MultiVector& porosity, 
      const Epetra_MultiVector& saturation, const Epetra_MultiVector& mol_density);

  void CalculateDiffusionTensor_(
      double md, int phase, const Epetra_MultiVector& porosity, 
      const Epetra_MultiVector& saturation, const Epetra_MultiVector& mol_density);

  int FindDiffusionValue(const std::string& tcc_name, double* md, int* phase);

  void CalculateAxiSymmetryDirection();

  // -- air-water partitioning using Henry's law. This is a temporary
  //    solution to get things moving.
  void PrepareAirWaterPartitioning_();
  void MakeAirWaterPartitioning_();

  // -- multiscale methods
  void AddMultiscalePorosity_(double t_old, double t_new, double t_int1, double t_int2);

  // initialization methods
  void InitializeAll_();
  void InitializeFieldFromField_(const std::string& field0, 
                                 const std::string& field1, 
                                 const Teuchos::Ptr<State>& S,
                                 bool call_evaluator, bool overwrite);

    // miscaleneous methods
    int FindComponentNumber(const std::string component_name);

    void ComputeVolumeDarcyFlux(Teuchos::RCP<const Epetra_MultiVector> flux,
            Teuchos::RCP<const Epetra_MultiVector> mol_den,
                                Teuchos::RCP<Epetra_MultiVector>& vol_darcy_flux);



 public:
    Teuchos::RCP<Teuchos::ParameterList> tp_list_;

    int MyPID;  // parallel information: will be moved to private
    int spatial_disc_order, temporal_disc_order, limiter_model;

    int nsubcycles;  // output information
    int internal_tests;
    double tests_tolerance;

 protected:
    Teuchos::RCP<TreeVector> soln_;

    Key domain_name_;
    Key saturation_key_;
    Key prev_saturation_key_;
    Key flux_key_;
    Key darcy_flux_key_;
    Key permeability_key_;
    Key tcc_key_;
    Key porosity_key_;
    Key tcc_matrix_key_;
    Key molar_density_key_;
    Key solid_residue_mass_key_;
    Key water_content_key_;
    Key mass_src_key_;
    Key conserve_qty_key_;
 
 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<State> S_;
  std::string passwd_;

  bool subcycling_, water_source_in_meters_;
  int dim;
  int saturation_name_;
  bool vol_flux_conversion_;

  Teuchos::RCP<CompositeVector> tcc_w_src;
  Teuchos::RCP<CompositeVector> tcc_tmp;  // next tcc
  Teuchos::RCP<CompositeVector> tcc;  // smart mirrow of tcc 
  Teuchos::RCP<Epetra_MultiVector> conserve_qty_, solid_qty_, water_qty_;
  Teuchos::RCP<const Epetra_MultiVector> flux_;
  Teuchos::RCP<const Epetra_MultiVector> ws_, ws_prev_, phi_, mol_dens_, mol_dens_prev_;
  Teuchos::RCP<Epetra_MultiVector> flux_copy_;
    
#ifdef ALQUIMIA_ENABLED
  Teuchos::RCP<AmanziChemistry::Alquimia_PK> chem_pk_;
  Teuchos::RCP<AmanziChemistry::ChemistryEngine> chem_engine_;
#endif

  Teuchos::RCP<Epetra_IntVector> upwind_cell_;
  Teuchos::RCP<Epetra_IntVector> downwind_cell_;

  Teuchos::RCP<const Epetra_MultiVector> ws_start, ws_end;  // data for subcycling 
  Teuchos::RCP<const Epetra_MultiVector> mol_dens_start, mol_dens_end;  // data for subcycling 
  Teuchos::RCP<Epetra_MultiVector> ws_subcycle_start, ws_subcycle_end;
  Teuchos::RCP<Epetra_MultiVector> mol_dens_subcycle_start, mol_dens_subcycle_end;

  int current_component_;  // data for lifting
  Teuchos::RCP<Operators::ReconstructionCell> lifting_;
  Teuchos::RCP<Operators::LimiterCell> limiter_;

  std::vector<Teuchos::RCP<TransportDomainFunction> > srcs_;  // Source or sink for components
  std::vector<Teuchos::RCP<TransportDomainFunction> > bcs_;  // influx BC for components
  double bc_scaling;
  Teuchos::RCP<Epetra_Vector> Kxy;  // absolute permeability in plane xy

  Teuchos::RCP<Epetra_Import> cell_importer;  // parallel communicators
  Teuchos::RCP<Epetra_Import> face_importer;

  // mechanical dispersion and molecual diffusion
  Teuchos::RCP<MDMPartition> mdm_;
  std::vector<WhetStone::Tensor> D_;

  bool flag_dispersion_;
  std::vector<int> axi_symmetry_;  // axi-symmetry direction of permeability tensor

  std::vector<Teuchos::RCP<MaterialProperties> > mat_properties_;  // vector of materials
  std::vector<Teuchos::RCP<DiffusionPhase> > diffusion_phase_;   // vector of phases

  // Hosting temporarily Henry law 
  bool henry_law_;
  std::vector<double> kH_;
  std::vector<int> air_water_map_;

  // multiscale models
  bool multiscale_porosity_;
  Teuchos::RCP<MultiscaleTransportPorosityPartition> msp_;

  double cfl_, dt_, dt_debug_, t_physics_;  

  std::vector<double> mass_solutes_exact_, mass_solutes_source_;  // mass for all solutes
  std::vector<double> mass_solutes_bc_, mass_solutes_stepstart_;
  std::vector<std::string> runtime_solutes_;  // names of trached solutes
  std::vector<std::string> runtime_regions_;

  int ncells_owned, ncells_wghost;
  int nfaces_owned, nfaces_wghost;
  int nnodes_wghost;
 
  std::vector<std::string> component_names_;  // details of components
  std::vector<double> mol_masses_;
  int num_aqueous, num_gaseous;
  double water_tolerance_, max_tcc_;
  bool dissolution_;

  // io
  Utils::Units units_;

  Teuchos::RCP<VerboseObject> vo_;
  Teuchos::RCP<Debugger> db_;

  // Forbidden.
  Transport_ATS(const Transport_ATS&);
  Transport_ATS& operator=(const Transport_ATS&);

 private:
  // factory registration
  static RegisteredPKFactory<Transport_ATS> reg_;
};

}  // namespace Transport
}  // namespace Amanzi

#endif

