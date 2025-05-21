/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Ethan Coon (coonet@ornl.gov)
*/

/*
  Transport PK

*/


/*!

This PK solves for transport of chemical species in water.  It may optionally
be paired with a PK for reactions, typically in an operator-split form, to
solve reactive transport.

The advection-dispersion equation for component *i* may be written as:

.. math::
  \frac{\partial (\Theta \xi_i)}{\partial t} = - \boldsymbol{\nabla} \cdot (\boldsymbol{q} \xi_i)
  + \boldsymbol{\nabla} \cdot (\Theta \, (\boldsymbol{D^*}_l + \tau \boldsymbol{D}_i) \boldsymbol{\nabla} \xi_i) + Q_s,

The primary variable for this PK is :math:`\xi_i`, the **mole ratio** of a
specie :math:`i` in the aqueous phase, with units of [mol i / mol H2O].

This seems a bit odd to most people who are used to reactive transport codes
where the primary variable is total component concentration :math:`C`, in units of [mol i
/ L].  The reason for this differences is to better enable variable-density
problems.  Note that the two are related:

.. math::
   C_i = n_l \xi_i

for molar density of liquid water :math:`n_l`.

For reactive transport problems, both concentration and mole ratio are output
in the visualization file.  For transport-only problems, concentration may be
output by adding a total component concentration evaluator that multiplies the
two quantities using an `"evaluator type`" = `"multiplicative evaluator`".

`"PK type`" = `"transport ATS`"

.. _pk-transport-ats-spec:
.. admonition:: pk-transport-ats-spec

   * `"domain name`" ``[string]`` **"domain"** specifies mesh name that defines
     the domain of this PK.

   * `"component names`" ``[Array(string)]`` **optional** No default. Provides
     the names of the components that will be transported.  Note that these
     must be provided by the user if transport is used alone, or are provided
     by the geochemical engine if reactions are to be used.

   * `"number of aqueous components`" ``[int]`` **-1** The total number of
     aqueous components.  Default value is the length of `"component names`"

   * `"boundary conditions`" ``[transport-bc-spec]`` Boundary conditions for
     transport are dependent on the boundary conditions for flow. See
     :ref:`Transport-specific Boundary Conditions`

   * `"molecular diffusivity [m^2 s^-1]`" ``[molecular-diffusivity-spec]`` See
     below.

   * `"tortuosity [-]`" ``[tortuosity-spec]`` See below.

   * `"source terms`" ``[transport-source-spec-list]`` Provides solute source.

   Math and solver algorithm options:

   * `"diffusion`" ``[pde-diffusion-typedinline-spec]`` Diffusion drives the
     distribution.  Typically we use finite volume here, but mimetic schemes
     may be used.  See :ref:`Diffusion`

   * `"diffusion preconditioner`" ``[pde-diffusion-typedinline-spec]`` Inverse
     of the above.  Likely only Jacobian term options are needed here, as the
     others default to the same as the `"diffusion`" list.  See
     :ref:`Diffusion`.

   * `"inverse`" ``[inverse-typed-spec]`` :ref:`Inverse` method for the
     diffusion-dispersion solve.  See :ref:`Inverse`.

   * `"cfl`" ``[double]`` **1.0** Courant-Friedrichs-Lewy condition, a limiter
     on the timestep size relative to spatial size.  Must be <= 1.

   * `"advection spatial discretization order`" ``[int]`` **1** Defines
     accuracy of the spatial discretization in the advection term.  It permits
     values 1 or 2. Default value is 1 (donor upwind) but 2 (a limiter scheme)
     is much less numerically diffusive, and recommended for most cases.

   * `"temporal discretization order`" ``[int]`` **1** Defines accuracy of
     temporal discretization.  It permits values 1 (forward Euler) and 2 (a
     Runga-Kutta scheme).

   * `"reconstruction`" ``[reconstruction-spec]`` collects reconstruction
     parameters for use in reconstructing the velocity field for 2nd order
     advection schemes.  See :ref:`Reconstructions`.

   KEYS

   - `"primary variable`" **"mole_fraction"** [mol C mol H2O^-1]
   - `"liquid water content`" **"water_content"** This variable is a multiplier
     in in the accumulation term. This is often just `"water_content`", but
     not in frozen cases, where it must be defined by the user (typically as
     :math:`\phi s n_l |V|` in the subsurface, or :math:`h n_l |V|` on the
     surface.
   - `"molar density liquid`" [mol H2O m^-3]
   - `"water flux`" The face-based water flux in [mol H2O s^-1].
   - `"water source`" Defines the water injection rate [mol H2O m^-2 s^-1] in
     surface and [mol H2O m^-3 s^-1] in subsurface) which applies to
     concentrations specified by the `"geochemical conditions`".  Note that if
     this PK is coupled to a surface flow PK, the unit of the water source
     there *must* be in [mol m^-2 s^-1], *not* in [m s^-1] as is an option for
     that PK (e.g. `"water source in meters`" must be set to `"false`" in the
     overland flow PK).
      
     The injection rate of a solute [molC s^-1], when given as the product of a
     concentration and a water source, is evaluated as:
      
        Concentration [mol C L^-1] :math:`\times`
        1000 [L m^-3] of water :math:`\times`
        water source [mol H2O m^-3 s^-1] :math:`\times`
        volume of injection domain [m^3] :math:`/`
        molar density of water [mol H2O m^-3]


Note, this is not dispersion, but strictly isotropic diffusion.        

.. _molecular-diffusivity-spec:
.. admonition:: molecular-diffusivity-spec

   For each aqueous component, a single value is provided for molecular diffusivity, e.g.

   `"COMPONENT_NAME`" ``[double]`` value [m^2 s^-1]


.. _tortuosity-spec:
.. admonition:: tortuosity-spec

   * `"aqueous`" ``[double]`` Defines tortuosity for calculating
     diffusivity of liquid solutes, [-].

   * `"gaseous`" ``[double]`` Defines tortuosity for calculating
     diffusivity of gas solutes, [-].  Not currently implemented!


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
#include "DenseVector.hh"
#include "PK_Physical_Default.hh"

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

namespace Amanzi {

// forward declarations
struct TensorVector;
namespace Operators {
class PDE_Diffusion;
class PDE_Accumulation;
class BCs;
class Operator;
}

namespace Transport {

// ummm -- why does this not use TreeVector? --ETC
class Transport_ATS : public PK_Physical_Default {
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

  void parseParameterList() override;

  // members required by PK interface
  virtual void Setup() override;

  // coupling with chemistry
#ifdef ALQUIMIA_ENABLED
  void setChemEngine(Teuchos::RCP<AmanziChemistry::Alquimia_PK> chem_pk);
#endif

  virtual void Initialize() override;

  virtual double get_dt() override;
  virtual void set_dt(double dt) override{};

  virtual bool AdvanceStep(double t_old, double t_new, bool reinit = false) override;
  virtual void CommitStep(double t_old, double t_new, const Tag& tag) override;
  virtual void CalculateDiagnostics(const Tag& tag) override{};

  // -- helper functions
  int get_num_aqueous_component() const {
    return num_aqueous_;
  }

 protected:

  // component-based or phase-based parameters
  using ParameterMap = std::map<std::string, double>;
  ParameterMap readParameterMapByComponent(Teuchos::ParameterList& plist, double default_val);
  ParameterMap readParameterMapByPhase(Teuchos::ParameterList& plist, double default_val);

  // transport physics members
  // -- calculation of a stable timestep needs saturations and darcy flux
  double ComputeStableTimeStep_();

  // -- setup/initialize helper functions
  void InitializeFields_();

  void SetupTransport_();
  void SetupPhysicalEvaluators_();

  // -- advance members -- portions of the operator
  void AdvanceAdvectionSources_RK1_(double t_old, double t_new, int spatial_order);
  void AdvanceAdvectionSources_RK2_(double t_old, double t_new, int spatial_order);
  void AdvanceDispersionDiffusion_(double t_old, double t_new);

  // -- helper functions used in AdvanceDispersionDiffusion_
  bool PopulateBoundaryData_(int component, Operators::BCs& bcs);

  // -- helper functions used in AdvanceAdvectionSources*
  void IdentifyUpwindCells_();
  void AddAdvection_FirstOrderUpwind_(double t_old,
          double t_new,
          const Epetra_MultiVector& tcc,
          Epetra_MultiVector& f);
  void AddAdvection_SecondOrderUpwind_(double t_old,
          double t_new,
          const Epetra_MultiVector& tcc,
          Epetra_MultiVector& f);
  void AddAdvection_SecondOrderUpwind_(double t_old,
          double t_new,
          const Epetra_Vector& tcc,
          Epetra_Vector& f,
          int component);

  void AddSourceTerms_(double t_old,
                       double t_new,
                       Epetra_MultiVector& conserve_qty,
                       int n0,
                       int n1);

  void InvertTccNew_(const Epetra_MultiVector& conserve_qty,
                     Epetra_MultiVector& tcc,
                     Epetra_MultiVector* solid_qty,
                     bool include_current_water_mass);

  // -- air-water partitioning using Henry's law. This is a temporary
  //    solution to get things moving.
  // void PrepareAirWaterPartitioning_();
  // void MakeAirWaterPartitioning_();

  // -- multiscale methods
  // void AddMultiscalePorosity_(double t_old, double t_new, double t_int1, double t_int2);

  // miscillaneous methods
  int FindComponentNumber_(const std::string& component_name) {
    auto comp = std::find(component_names_.begin(), component_names_.end(), component_name);
    if (comp == component_names_.end()) return -1;
    else return comp - component_names_.begin();
  }

 protected:

  // dependencies
  Key flux_key_; // advecting water flux [mol / s] on faces
  Key lwc_key_; // liquid water content [mol]
  Key cv_key_;

  // workspace
  Key solid_residue_mass_key_; // residue -- mass that was left behind by water
                               // sinks that do not carry aqueous species, e.g. freezing or evaporation
  Key conserve_qty_key_;
  Key dispersion_tensor_key_;

  // component information
  std::vector<std::string> component_names_; // details of components
  int num_components_, num_aqueous_;

  // parameters
  ParameterMap molar_masses_; // molar mass [kg / mol C] by component
  ParameterMap tcc_max_; // maximum valid concentration [mol C / mol H2O] by component

  ParameterMap molec_diff_; // molecular diffusivity by component [m^2 s^-1]
  ParameterMap tortuosity_; // phase-based tortuosity [-]
  double water_tolerance_; // mol H2O / m^d (d = 2 for surface, d = 3 for subsurface)

#ifdef ALQUIMIA_ENABLED
  Teuchos::RCP<AmanziChemistry::Alquimia_PK> chem_pk_;
  Teuchos::RCP<AmanziChemistry::ChemistryEngine> chem_engine_;
#else
  Teuchos::RCP<bool> chem_pk_;
  Teuchos::RCP<bool> chem_engine_;
#endif

  // temporal integration
  // -- discretization order: RK1 (forward Euler) or RK2 (predictor-corrector)
  int temporal_disc_order_;
  double cfl_;
  double dt_stable_, dt_max_;

  // Source terms
  // srcs should go away, and instead use vectors from State and evaluators --ETC
  Key water_src_key_;
  Key geochem_src_factor_key_;
  bool has_water_src_key_;
  bool water_src_in_meters_;
  Key molar_dens_key_;
  std::vector<Teuchos::RCP<TransportDomainFunction>> srcs_; // Source or sink for components

  // operators for advection
  int adv_spatial_disc_order_;
  Teuchos::RCP<Epetra_IntVector> upwind_cell_, downwind_cell_;
  Teuchos::RCP<Operators::ReconstructionCellLinear> lifting_;
  Teuchos::RCP<Operators::LimiterCell> limiter_;

  Teuchos::RCP<Operators::BCs> adv_bcs_;
  std::vector<Teuchos::RCP<TransportDomainFunction>> bcs_;

  // operators for dispersion/diffusion
  bool has_diffusion_, has_dispersion_;
  Teuchos::RCP<TensorVector> D_; // workspace, disp + diff
  Teuchos::RCP<Operators::BCs> diff_bcs_;
  Teuchos::RCP<Operators::Operator> diff_global_op_;
  Teuchos::RCP<Operators::PDE_Diffusion> diff_op_;
  Teuchos::RCP<Operators::PDE_Accumulation> diff_acc_op_;
  Teuchos::RCP<CompositeVector> diff_sol_; // workspace

  // io
  Utils::Units units_;

 private:
  // Forbidden.
  Transport_ATS(const Transport_ATS&) = delete;
  Transport_ATS& operator=(const Transport_ATS&) = delete;

 private:
  // factory registration
  static RegisteredPKFactory<Transport_ATS> reg_;
};

// helper functions
void CheckGEDProperty(const Epetra_MultiVector& tracer,
                      double t_physics);

void CheckTracerBounds(const Epetra_MultiVector& tcc,
                       const Epetra_MultiVector& tcc_prev,
                       const AmanziMesh::Mesh& mesh,
                       double t_physics,
                       int component,
                       double lower_bound,
                       double upper_bound,
                       double tol = 0.0);

} // namespace Transport
} // namespace Amanzi

#endif
