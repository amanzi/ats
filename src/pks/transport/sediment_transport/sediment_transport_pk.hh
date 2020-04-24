/*
  Transport PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Daniil Svyatsky (dasvyat@lanl.gov)
*/

#ifndef AMANZI_ATS_SEDIMENTTRANSPORT_PK_HH_
#define AMANZI_ATS_SEDIMENTTRANSPORT_PK_HH_

// TPLs
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Import.h"
#include "Teuchos_RCP.hpp"

// Amanzi
#include "CompositeVector.hh"
#include "Explicit_TI_FnBase.hh"
//#include "MaterialProperties.hh"
#include "PK.hh"
#include "PK_Factory.hh"
#include "ReconstructionCell.hh"
#include "State.hh"
#include "Tensor.hh"
#include "Units.hh"
#include "VerboseObject.hh"
#include "PK_PhysicalExplicit.hh"
#include "DenseVector.hh"

#include <string>

// Transport
#include "TransportDomainFunction.hh"
#include "SedimentTransportDefs.hh"


/* ******************************************************************
The transport PK receives a reduced (optional) copy of a physical 
state at time n and returns a different state at time n+1. 

Unmodified physical quantaties in the returned state are the smart 
pointers to the original variables.
****************************************************************** */

namespace Amanzi {
namespace SedimentTransport {

typedef double AnalyticFunction(const AmanziGeometry::Point&, const double);

  class SedimentTransport_PK : public PK_PhysicalExplicit<Epetra_Vector> {

  public:
    SedimentTransport_PK(Teuchos::ParameterList& pk_tree,
               const Teuchos::RCP<Teuchos::ParameterList>& glist,
               const Teuchos::RCP<State>& S,
               const Teuchos::RCP<TreeVector>& soln);

    SedimentTransport_PK(const Teuchos::RCP<Teuchos::ParameterList>& glist,
                 Teuchos::RCP<State> S,
                 const std::string& pk_list_name,
                 std::vector<std::string>& component_names);

    ~SedimentTransport_PK() = default;

  // members required by PK interface
  virtual void Setup(const Teuchos::Ptr<State>& S);
  virtual void Initialize(const Teuchos::Ptr<State>& S);

  virtual double get_dt();
  virtual void set_dt(double dt) {};

  virtual bool AdvanceStep(double t_old, double t_new, bool reinit=false); 
  virtual void CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S);
  virtual void CalculateDiagnostics(const Teuchos::RCP<State>& S) {};

  virtual void set_states(const Teuchos::RCP<const State>& S,
                          const Teuchos::RCP<State>& S_inter,
                          const Teuchos::RCP<State>& S_next);

  virtual std::string name() { return "sediment transport"; }
  Key get_domain_name() {return domain_name_;}

  // main transport members
  // -- calculation of a stable time step needs saturations and darcy flux
  double StableTimeStep();
  void Sinks2TotalOutFlux(Epetra_MultiVector& tcc,
                          std::vector<double>& total_outflux, int n0, int n1);

  // -- access members  
  inline double cfl() { return cfl_; }
  Teuchos::RCP<const State> state() { return S_; }
  Teuchos::RCP<CompositeVector> total_component_concentration() { return tcc_tmp; }

  // -- control members
  void CreateDefaultState(Teuchos::RCP<const AmanziMesh::Mesh>& mesh, int ncomponents);
  // void Policy(Teuchos::Ptr<State> S);

  // void VV_CheckGEDproperty(Epetra_MultiVector& tracer) const; 
  // void VV_CheckTracerBounds(Epetra_MultiVector& tracer, int component,
  //                           double lower_bound, double upper_bound, double tol = 0.0) const;
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

  bool PopulateBoundaryData(std::vector<int>& bc_model,
                            std::vector<double>& bc_value, int component);

  // -- limiters 
  void LimiterBarthJespersen(const int component,
                             Teuchos::RCP<const Epetra_Vector> scalar_field, 
                             Teuchos::RCP<CompositeVector>& gradient, 
                             Teuchos::RCP<Epetra_Vector>& limiter);

  // const std::vector<std::string>  component_names(){return component_names_;};
  // int num_aqueous_component() {return num_aqueous;};


 private:
  void InitializeFields_(const Teuchos::Ptr<State>& S);

  // advection members
  void AdvanceDonorUpwind(double dT);
  // void AdvanceSecondOrderUpwindRKn(double dT);
  // void AdvanceSecondOrderUpwindRK1(double dT);
  // void AdvanceSecondOrderUpwindRK2(double dT);
  void Advance_Diffusion(double t_old, double t_new);

  // time integration members
    void FunctionalTimeDerivative(const double t, const Epetra_Vector& component, Epetra_Vector& f_component){};
    //  void Functional(const double t, const Epetra_Vector& component, TreeVector& f_component);

  void IdentifyUpwindCells();

  void InterpolateCellVector(
      const Epetra_MultiVector& v0, const Epetra_MultiVector& v1, 
      double dT_int, double dT, Epetra_MultiVector& v_int);

  const Teuchos::RCP<Epetra_IntVector>& upwind_cell() { return upwind_cell_; }
  const Teuchos::RCP<Epetra_IntVector>& downwind_cell() { return downwind_cell_; }  

  // physical models
  // -- dispersion and diffusion
  void CalculateDiffusionTensor_(const Epetra_MultiVector& km,
                                 const Epetra_MultiVector& saturation,
                                 const Epetra_MultiVector& mol_density);

  int FindDiffusionValue(const std::string& tcc_name, double* md, int* phase);

  void CalculateAxiSymmetryDirection();

  // initialization methods
  void InitializeAll_();
  void InitializeFieldFromField_(const std::string& field0, 
                                 const std::string& field1, 
                                 const Teuchos::Ptr<State>& S,
                                 bool call_evaluator, bool overwrite);

 public:
    Teuchos::RCP<Teuchos::ParameterList> tp_list_;
    Teuchos::RCP<const Teuchos::ParameterList> preconditioner_list_;
    Teuchos::RCP<const Teuchos::ParameterList> linear_solver_list_;
    Teuchos::RCP<const Teuchos::ParameterList> nonlinear_solver_list_;

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
    Key tcc_key_;
    Key molar_density_key_;
    Key solid_residue_mass_key_;
    Key sd_trapping_key_, sd_settling_key_, sd_erosion_key_, horiz_mixing_key_, porosity_key_, sd_organic_key_;
    Key elevation_increase_key_;

  
 
 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<State> S_;
  std::string passwd_;

  bool subcycling_;
  int dim;
  int saturation_name_;
  bool vol_flux_conversion_;

  Teuchos::RCP<CompositeVector> tcc_tmp;  // next tcc
  Teuchos::RCP<CompositeVector> tcc;  // smart mirrow of tcc 
  Teuchos::RCP<Epetra_MultiVector> conserve_qty_, solid_qty_;
  Teuchos::RCP<const Epetra_MultiVector> flux_;
  Teuchos::RCP<const Epetra_MultiVector> ws_, ws_prev_, mol_dens_;//, mol_dens_prev_;
  Teuchos::RCP<Epetra_MultiVector> flux_copy_;
  Teuchos::RCP<const Epetra_MultiVector> km_;  
    
  Teuchos::RCP<Epetra_IntVector> upwind_cell_;
  Teuchos::RCP<Epetra_IntVector> downwind_cell_;

  Teuchos::RCP<const Epetra_MultiVector> ws_start, ws_end;  // data for subcycling 
  Teuchos::RCP<const Epetra_MultiVector> mol_dens_start, mol_dens_end;  // data for subcycling 
  Teuchos::RCP<Epetra_MultiVector> ws_subcycle_start, ws_subcycle_end;
  Teuchos::RCP<Epetra_MultiVector> mol_dens_subcycle_start, mol_dens_subcycle_end;

  int current_component_;  // data for lifting
  Teuchos::RCP<Operators::ReconstructionCell> lifting_;

  std::vector<Teuchos::RCP<TransportDomainFunction> > srcs_;  // Source or sink for components
  std::vector<Teuchos::RCP<TransportDomainFunction> > bcs_;  // influx BC for components
  double bc_scaling;

  Teuchos::RCP<Epetra_Import> cell_importer;  // parallel communicators
  Teuchos::RCP<Epetra_Import> face_importer;

  // mechanical dispersion and molecual diffusion
  // Teuchos::RCP<MDMPartition> mdm_;
    
  std::vector<WhetStone::Tensor> D_;
  std::string diffusion_preconditioner, diffusion_solver;    

  // bool flag_dispersion_;
  // std::vector<int> axi_symmetry_;  // axi-symmetry direction of permeability tensor
  

  // std::vector<Teuchos::RCP<MaterialProperties> > mat_properties_;  // vector of materials
  // std::vector<Teuchos::RCP<DiffusionPhase> > diffusion_phase_;   // vector of phases


  double cfl_, dt_, dt_debug_, t_physics_;  

  double mass_sediment_exact_, mass_sediment_source_;  // mass for all sediment
  double mass_sediment_bc_, mass_sediment_stepstart_;
  std::vector<std::string> runtime_sediment_;  // names of trached sediment
  std::vector<std::string> runtime_regions_;

  int ncells_owned, ncells_wghost;
  int nfaces_owned, nfaces_wghost;
  int nnodes_wghost;
 
  std::vector<double> mol_masses_;
  double sediment_density_;
  int num_aqueous;

  std::vector<std::string> component_names_;  // details of components   
  double water_tolerance_, max_tcc_;  

  // io
    Utils::Units units_;
    Teuchos::RCP<VerboseObject> vo_;

  // Forbidden.
  SedimentTransport_PK(const SedimentTransport_PK&);
  SedimentTransport_PK& operator=(const SedimentTransport_PK&);

 private:
  // factory registration
  static RegisteredPKFactory<SedimentTransport_PK> reg_;
};

}  // namespace SedimentTransport
}  // namespace Amanzi

#endif

