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


/*!

The advection-diffusion equation for component *i* in partially saturated porous media may be written as

.. math::
  \frac{\partial (\phi s_l C_i)}{\partial t}
  =
  - \boldsymbol{\nabla} \cdot (\boldsymbol{q} C_i)
  + \boldsymbol{\nabla} \cdot (\phi s_l\, (\boldsymbol{D^*}_l + \tau \boldsymbol{D}_i) \boldsymbol{\nabla} C_i) + Q_s,

The advection-diffusion equation for component *i* in the surface may be written as

.. math::
  \frac{\partial (C_i)}{\partial t}
  =
  - \boldsymbol{\nabla} \cdot (\boldsymbol{q_s} C_i)
  + \boldsymbol{\nabla} \cdot ( (\boldsymbol{D^*}_l + \tau \boldsymbol{D}_i) \boldsymbol{\nabla} C_i) + Q_s,

.. _transport-spec:
.. admonition:: transport-spec

    * `"PK type`" ``[string]`` **"transport ats"**

    * `"domain name`" ``[string]`` **domain** specifies mesh name that defines
      the domain of this PK.

    * `"component names`" ``[Array(string)]`` No default. Provides the names of the
      components that will be transported. Must be in the order: aqueous, gaseous, solid.

    * `"number of aqueous components`" ``[int]`` **-1** The total number of
      aqueous components.  Default value is the length of `"component names`"

    * `"number of gaseous components`" ``[int]`` **0** The total number of
      gaseous components.

    * `"boundary conditions`" ``[transport-bc-spec]`` Boundary conditions for
      transport are dependent on the boundary conditions for flow. See
      `Flow-specific Boundary Conditions`_ and `Transport-specific Boundary Conditions`_

    * `"component molar masses`" ``[Array(double)]`` No default. Molar mass of
      each component.

    * `"molecular diffusion`" ``[molecular-diffusion-spec]`` defines names of
      solutes in aqueous and gaseous phases and related diffusivity values.

    * "material properties" [material-properties-spec-list] Defines material
      properties see below).

    Source terms:

    * `"source terms`" [transport-source-spec-list] Provides solute source.

    Physical model and assumptions:

    * `"physical models and assumptions`" [material-properties-spec] Defines material properties.

    * `"effective transport porosity`" [bool] If *true*, effective transport porosity
      will be used by dispersive-diffusive fluxes instead of total porosity.
      Default is *false*.

    Math and solver algorithm options:

    * `"diffusion`" ``[pde-diffusion-spec]`` Diffusion drives the distribution.
      Typically we use finite volume here.  See PDE_Diffusion_

    * `"diffusion preconditioner`" ``[pde-diffusion-spec]`` Inverse of the
      above.  Likely only Jacobian term options are needed here, as the others
      default to the same as the `"diffusion`" list.  See PDE_Diffusion_.

    * `"inverse`" ``[inverse-typed-spec]`` Inverse_ method for the solve.

    * `"cfl`" [double] Time step limiter, a number less than 1. Default value is 1.

    * `"spatial discretization order`" [int] defines accuracy of spatial discretization.
      It permits values 1 or 2. Default value is 1.

    * `"temporal discretization order`" [int] defines accuracy of temporal discretization.
      It permits values 1 or 2 and values 3 or 4 when expert parameter
      `"generic RK implementation`" is set to true. Note that RK3 is not monotone.
      Default value is 1.

    * `"reconstruction`" [list] collects reconstruction parameters. The available options are
      describe in the separate section below.

    * `"transport subcycling`" ``[bool]`` **true** The code will default to
      subcycling for transport within the master PK if there is one.


    Developer parameters:

    * `"enable internal tests`" [bool] turns on various internal tests during
      run time. Default value is `"false`".

    * `"generic RK implementation`" [bool] leads to generic implementation of
      all Runge-Kutta methods. Default value is `"false`".

    * `"internal tests tolerance`" [double] tolerance for internal tests such as the
      divergence-free condition. The default value is 1e-6.

    * `"runtime diagnostics: solute names`" [Array(string)] defines solutes that will be
      tracked closely each time step if verbosity `"high`". Default value is the first
      solute in the global list of `"aqueous names`" and the first gas in the global list
      of `"gaseous names`".

    * `"runtime diagnostics: regions`" [Array(string)] defines a boundary region for
      tracking solutes. Default value is a seepage face boundary, see Flow PK.

    KEYS

    - `"saturation liquid`" This variable is a multiplier in in the
      accumulation term. For subsurface transport, this will typically be the
      saturation (`"saturation_liquid`"). For surface transport, this will
      typically be the ponded depth (`"ponded_depth`").

    - `"previous saturation liquid`"

    - `"molar density liquid`"  Transport is solved
      for concentrations in units of mol fractions. Molar density is needed for conversion.

    - `"water flux`"

    - `"water source`" Defines the water injection rate [mol H2O m^-2 s^-1] in
      surface and [mol H2O m^-3 s^-1] in subsurface) which applies to
      concentrations specified by the `"geochemical conditions`".  Note that if
      this PK is coupled to a surface flow PK, the unit of the water source
      there *must* be in [mol m^-2 s^-1], *not* in [m s^-1] as is an option for
      that PK (e.g. `"water source in meters`" must be set to `"false`" in the
      overland flow PK).

      The injection rate of a solute [molC s^-1], when given as the product of
      a concentration and a water source, is evaluated as:

      Concentration [mol C L^-1] *
        1000 [L m^-3] of water *
        water source [mol H2O m^-3 s^-1] *
        volume of injection domain [m^3] /
        molar density of water [mol H2O m^-3]


.. _molecular-diffusion-spec:
.. admonition:: molecular-diffusion-spec

   * `"aqueous names`" ``[Array(string)]`` List of aqueous component names to
     be diffused.
   * `"aqueous values`" ``[Array(string)]`` Diffusivities of each component.


.. code-block:: xml

   <ParameterList name="molecular diffusion">
     <Parameter name="aqueous names" type=Array(string)" value="{CO2(l),Tc-99}"/>
     <Parameter name="aqueous values" type=Array(double)" value="{1e-8,1e-9}"/>
   </ParameterList>


.. _material-properties-spec:
.. admonition:: material-properties-spec

   * `"region`" ``[Array(string)]`` Defines geometric regions for material SOIL.

   * `"model`" ``[string]`` **scalar** Defines dispersivity model.  One of:

     - `"scalar`" : scalar dispersivity
     - `"Bear`" : dispersion split into along- and across- velocity
     - `"Burnett-Frind`"
     - `"Lichtner-Kelkar-Robinson`"

   * `"parameters for MODEL`" ``[list]`` where `"MODEL`" is the model name.

   IF model == scalar

   ONE OF

   * `"alpha`" ``[double]`` defines dispersivity in all directions, [m].

   OR

   * `"dispersion coefficient`" ``[double]`` defines dispersion coefficient [m^2 s^-1].

   END

   ELSE IF model == Bear

   * `"alpha_l`" ``[double]`` defines dispersion in the direction of Darcy velocity, [m].
   * `"alpha_t`" ``[double]`` defines dispersion in the orthogonal direction, [m].

   ELSE IF model == Burnett-Frind

   * `"alphaL`" ``[double]`` defines the longitudinal dispersion in the direction
     of Darcy velocity, [m].
   * `"alpha_th`" ``[double]`` Defines the transverse dispersion in the horizonal
     direction orthogonal directions, [m].
   * `"alpha_tv`" ``[double]`` Defines dispersion in the orthogonal directions,
     [m].  When `"alpha_th`" equals to `"alpha_tv`", we obtain dispersion in
     the direction of the Darcy velocity.

   ELSE IF model == Lichtner-Kelker-Robinson

   * `"alpha_lh`" ``[double]`` defines the longitudinal dispersion in the
     horizontal direction, [m].
   * `"alpha_lv`" ``[double]`` Defines the longitudinal dispersion in the vertical
     direction, [m].  When `"alpha_lh`" equals to `"alpha_lv`", we obtain
     dispersion in the direction of the Darcy velocity.
   * `"alpha_th`" ``[double]`` Defines the transverse dispersion in the horizontal
     direction orthogonal directions, [m].
   * `"alpha_tv" ``[double]`` Defines dispersion in the orthogonal directions.
     When `"alpha_th`" equals to `"alpha_tv`", we obtain dispersion in the
     direction of the Darcy velocity.

   END

   * `"aqueous tortuosity`" ``[double]`` Defines tortuosity for calculating
     diffusivity of liquid solutes, [-].

   * `"gaseous tortuosity`" ``[double]`` Defines tortuosity for calculating
     diffusivity of gas solutes, [-].


.. _transport-source-spec:
.. admonition:: transport-source-spec

   * `"component mass source`" ``[list]``  Defines solute source injection rate.

     * `"spatial distribution method`" ``[string]`` One of:

        - `"volume`", source is considered as extensive quantity [molC s^-1] and is evenly distributed across the region.
        - `"none`", source is considered as intensive quantity. [molC m^-2 s^-1] in surface and [molC m^-3 s^-1] in subsurface

     * `"geochemical`" ``[list]``  Defines a source by setting solute concentration for all components (in moles/L) and an injection
       rate given by the water source.  Currently, this option is only available for Alquimia provided geochemical conditions.

       - `"geochemical conditions`" ``[Array(string)]`` List of geochemical constraints providing concentration for solute injection.

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
#include "ReconstructionCellLinear.hh"
#include "State.hh"
#include "Tensor.hh"
#include "Units.hh"
#include "VerboseObject.hh"
#include "Debugger.hh"
#include "PK_PhysicalExplicit.hh"
#include "DenseVector.hh"

#include <string>

#ifdef ALQUIMIA_ENABLED
#  include "Alquimia_PK.hh"
#  include "ChemistryEngine.hh"
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

typedef double
AnalyticFunction(const AmanziGeometry::Point&, const double);

// ummm -- why does this not use TreeVector? --ETC
class Transport_ATS : public PK_PhysicalExplicit<Epetra_Vector> {
 public:
  Transport_ATS(Teuchos::ParameterList& pk_tree,
                const Teuchos::RCP<Teuchos::ParameterList>& glist,
                const Teuchos::RCP<State>& S,
                const Teuchos::RCP<TreeVector>& solution);

  Transport_ATS(const Teuchos::RCP<Teuchos::ParameterList>& glist,
                Teuchos::RCP<State> S,
                const std::string& pk_list_name,
                std::vector<std::string>& component_names);

  ~Transport_ATS() = default;

  // members required by PK interface
  virtual void Setup() override;
  virtual void Initialize() override;

  virtual double get_dt() override;
  virtual void set_dt(double dt) override{};
  virtual void set_tags(const Tag& current, const Tag& next) override;

  virtual bool AdvanceStep(double t_old, double t_new, bool reinit = false) override;
  virtual void CommitStep(double t_old, double t_new, const Tag& tag) override;
  virtual void CalculateDiagnostics(const Tag& tag) override{};

  // main transport members
  // -- calculation of a stable time step needs saturations and darcy flux
  double StableTimeStep();
  void
  Sinks2TotalOutFlux(Epetra_MultiVector& tcc, std::vector<double>& total_outflux, int n0, int n1);

  // coupling with chemistry
#ifdef ALQUIMIA_ENABLED
  void SetupAlquimia(Teuchos::RCP<AmanziChemistry::Alquimia_PK> chem_pk,
                     Teuchos::RCP<AmanziChemistry::ChemistryEngine> chem_engine);
#endif

  // -- access members
  inline double get_cfl() { return cfl_; }
  Teuchos::RCP<const State> get_state() { return S_; }
  Teuchos::RCP<CompositeVector> get_total_component_concentration() { return tcc_tmp; }

  // -- control members
  void CreateDefaultState(Teuchos::RCP<const AmanziMesh::Mesh>& mesh, int ncomponents);
  void Policy(const Tag& tag);

  void VV_CheckGEDproperty(Epetra_MultiVector& tracer) const;
  void VV_CheckTracerBounds(Epetra_MultiVector& tracer,
                            int component,
                            double lower_bound,
                            double upper_bound,
                            double tol = 0.0) const;
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
  void ComputeAddSourceTerms(double tp, double dtp, Epetra_MultiVector& tcc, int n0, int n1);

  // void MixingSolutesWthSources(double told, double tnew);

  bool
  PopulateBoundaryData(std::vector<int>& bc_model, std::vector<double>& bc_value, int component);

  // -- limiters
  void LimiterBarthJespersen(const int component,
                             Teuchos::RCP<const Epetra_Vector> scalar_field,
                             Teuchos::RCP<CompositeVector>& gradient,
                             Teuchos::RCP<Epetra_Vector>& limiter);

  const std::vector<std::string> get_component_names() { return component_names_; };
  int get_num_aqueous_component() { return num_aqueous; };
  int get_num_gaseous_component() { return num_gaseous; };


 private:
  void InitializeFields_();

  // advection members
  void AdvanceDonorUpwind(double dT);
  void AdvanceSecondOrderUpwindRKn(double dT);
  void AdvanceSecondOrderUpwindRK1(double dT);
  void AdvanceSecondOrderUpwindRK2(double dT);
  void Advance_Dispersion_Diffusion(double t_old, double t_new);

  // time integration members
  void FunctionalTimeDerivative(const double t,
                                const Epetra_Vector& component,
                                Epetra_Vector& f_component) override;
  //  void FunctionalTimeDerivative(const double t, const Epetra_Vector& component, TreeVector& f_component);

  void IdentifyUpwindCells();

  void InterpolateCellVector(const Epetra_MultiVector& v0,
                             const Epetra_MultiVector& v1,
                             double dT_int,
                             double dT,
                             Epetra_MultiVector& v_int);

  const Teuchos::RCP<Epetra_IntVector>& get_upwind_cell() { return upwind_cell_; }
  const Teuchos::RCP<Epetra_IntVector>& get_downwind_cell() { return downwind_cell_; }

  // physical models
  // -- dispersion and diffusion
  void CalculateDispersionTensor_(const Epetra_MultiVector& darcy_flux,
                                  const Epetra_MultiVector& porosity,
                                  const Epetra_MultiVector& saturation,
                                  const Epetra_MultiVector& mol_density);

  void CalculateDiffusionTensor_(double md,
                                 int phase,
                                 const Epetra_MultiVector& porosity,
                                 const Epetra_MultiVector& saturation,
                                 const Epetra_MultiVector& mol_density);

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
  void InitializeFieldFromField_(const Key& field0,
                                 const Tag& tag0,
                                 const Key& field1,
                                 const Tag& tag1,
                                 bool call_evaluator,
                                 bool overwrite);

  // miscaleneous methods
  int FindComponentNumber(const std::string component_name);

  void ComputeVolumeDarcyFlux(Teuchos::RCP<const Epetra_MultiVector> flux,
                              Teuchos::RCP<const Epetra_MultiVector> mol_den,
                              Teuchos::RCP<Epetra_MultiVector>& vol_darcy_flux);

 public:
  int MyPID; // parallel information: will be moved to private
  int spatial_disc_order, temporal_disc_order, limiter_model;

  int nsubcycles; // output information
  int internal_tests;
  double tests_tolerance;

 protected:
  Key saturation_key_;
  Key flux_key_;
  Key darcy_flux_key_;
  Key permeability_key_;
  Key tcc_key_;
  Key porosity_key_;
  Key tcc_matrix_key_;
  Key molar_density_key_;
  Key solid_residue_mass_key_;
  Key water_content_key_;
  Key water_src_key_, solute_src_key_;
  bool has_water_src_key_;
  bool water_src_in_meters_;
  Key geochem_src_factor_key_;
  Key conserve_qty_key_;
  Key cv_key_;

 private:
  bool subcycling_;
  int dim;
  int saturation_name_;
  bool vol_flux_conversion_;

  Key passwd_;

  Teuchos::RCP<CompositeVector> tcc_w_src;
  Teuchos::RCP<CompositeVector> tcc_tmp; // next tcc
  Teuchos::RCP<CompositeVector> tcc;     // smart mirrow of tcc
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

  Teuchos::RCP<const Epetra_MultiVector> ws_current, ws_next;             // data for subcycling
  Teuchos::RCP<const Epetra_MultiVector> mol_dens_current, mol_dens_next; // data for subcycling
  Teuchos::RCP<Epetra_MultiVector> ws_subcycle_current, ws_subcycle_next;
  Teuchos::RCP<Epetra_MultiVector> mol_dens_subcycle_current, mol_dens_subcycle_next;

  int current_component_; // data for lifting
  Teuchos::RCP<Operators::ReconstructionCellLinear> lifting_;
  Teuchos::RCP<Operators::LimiterCell> limiter_;

  std::vector<Teuchos::RCP<TransportDomainFunction>> srcs_; // Source or sink for components
  std::vector<Teuchos::RCP<TransportDomainFunction>> bcs_;  // influx BC for components
  double bc_scaling;
  Teuchos::RCP<Epetra_Vector> Kxy; // absolute permeability in plane xy

  Teuchos::RCP<Epetra_Import> cell_importer; // parallel communicators
  Teuchos::RCP<Epetra_Import> face_importer;

  // mechanical dispersion and molecual diffusion
  Teuchos::RCP<MDMPartition> mdm_;
  std::vector<WhetStone::Tensor> D_;

  bool flag_dispersion_;
  std::vector<int> axi_symmetry_; // axi-symmetry direction of permeability tensor

  std::vector<Teuchos::RCP<MaterialProperties>> mat_properties_; // vector of materials
  std::vector<Teuchos::RCP<DiffusionPhase>> diffusion_phase_;    // vector of phases

  // Hosting temporarily Henry law
  bool henry_law_;
  std::vector<double> kH_;
  std::vector<int> air_water_map_;

  double cfl_, dt_, dt_debug_, t_physics_;

  std::vector<double> mass_solutes_exact_, mass_solutes_source_; // mass for all solutes
  std::vector<double> mass_solutes_bc_, mass_solutes_stepstart_;
  std::vector<std::string> runtime_solutes_; // solutes tracked for diagnostics
  std::vector<std::string> runtime_regions_;

  int ncells_owned, ncells_wghost;
  int nfaces_owned, nfaces_wghost;
  int nnodes_wghost;

  std::vector<std::string> component_names_; // details of components
  std::vector<double> mol_masses_;
  int num_aqueous, num_gaseous, num_components, num_primary, num_advect;
  double water_tolerance_, max_tcc_;
  bool dissolution_;

  // io
  Utils::Units units_;
  Tag tag_subcycle_;
  Tag tag_subcycle_current_;
  Tag tag_subcycle_next_;
  Tag tag_flux_next_ts_; // what is this? --ETC

 private:
  // Forbidden.
  Transport_ATS(const Transport_ATS&) = delete;
  Transport_ATS& operator=(const Transport_ATS&) = delete;

 private:
  // factory registration
  static RegisteredPKFactory<Transport_ATS> reg_;
};

} // namespace Transport
} // namespace Amanzi

#endif
