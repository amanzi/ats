/*--------------------------------------------------------------------------
  ATS

  License: see $ATS_DIR/COPYRIGHT
  Author: Andrew Graus (agraus@lbl.gov)

  --------------------------------------------------------------------------*/
/*!

This PK couples ATS to the BGC code EcoSIM. This PK essentially takes over the
surface balance aspects of ATS and replaces them with those from EcoSIM.

In addition this code takes data from ATS state and loads it into a struct that
is fortran readable (BGCContainer). The data structures and methods were adapted
from those used in the Alquimia code, additionally the general code structure
Engine code are adapted from Alquimia as well.

Structures for looping over cells of columns were adapted from ATS's simpleBGC code

`"PK type`" = `"EcoSIM for ATS`"

.. _pk-ecosim-spec:
.. admonition:: pk-ecosim-spec

   * `"engine`" ``[string]`` **EcoSIM**  Engine name - inspired by Alquimias options if this
     code is adapted to used to drive other BGC codes.

   * `"heat capacity [MJ mol^-1 K^-1]`" ``[double]]`` **0.02** heat capacity of the soil layers

   * `"Field Capacity [MPa]`" ``[double]]`` **-0.033** pressure at field capacity

   * `"Wilting Point [MPa]`" ``[double]]`` **-1.5** pressure at wilting point

   * `"initial time step [s]`" ``[double]]`` **3600.0** EcoSIM is an hourly model so the
     standard is to run it after the end of hour 1

   * `"EcoSIM Precipitation`" ``[bool]`` This allows EcoSIM to partition the precipitation
     itself. If false it will expect the precipitation forcing to be already divided into
     rain and snow as in ATS.

   * `"Prescribe Albedo`" ``[bool]`` determines if the code will use EcoSIM's internal
     albedo calculation, or if it will be prescirbed from data.

   * `"Prescribe Phenology`" ``[bool]`` Determines if the coupling will use the prescribed
     phenology methodology where LAI and PFT are input. This will eventually be complemented
     with a full phenology method, where plant parameters are directly, but will remain in
     the code as an option.

   * `"starting day of year [0-364]`" ``[int]`` day of the year, needed for EcoSIMs internal
     radiation and biogeochemical processes.

   * `"Starting year`" ``[int]`` Year also needed for internal EcoSIM computations

    * `"Number of PFTs [1-5] "`" ``[int]`` Number of PFTs allowed in every column. 5
     is the maximum number allowed

   * `"domain name`" ``[string]`` **domain**

   * `"surface domain name`" ``[string]`` **surface**

   EVALUATORS:

   - `"Bulk Density`" `[ ]`
   - `"Hydraulic Conductivity `[Pa]`
   - `"Matric Pressure`" `[Pa]`

   DEPENDENCIES
   //Sources
   `"surface water source ecosim`"     **surface-ecosim_water_source**
   `"surface energy source ecosim`"    **surface-ecosim_source**
   `"subsurface water source ecosim`"  **ecosim_water_source**
   `"surface water source ecosim`"     **surface-ecosim_water_source**

   //surface balance variables
   `"incoming shortwave radiation`"      **surface-incoming_shortwave_radiation**
   `"incoming longwave radiation`"       **surface-incoming_longwave_radiation**
   `"air temperature`"                   **surface-air_temperature**
   `"vapor pressure air`"                **surface-vapor_pressure_air**
   `"wind speed`"                        **surface-wind_speed**
   `"precipitation rain`"                **surface-precipitation_rain**
   `"precipitation snow`"                **surface-precipitation_snow**
   `"precipitation total`"               **surface-precipitation_total**
   `"snow depth`"                        **surface-snow_depth**
   `"snow_albedo`"                       **surface-snow_albedo**
   `"snow temperature`"                  **surface-snow_temperature**
   `"canopy longwave radiation`"         **surface-canopy_longwave_radiation**
   `"canopy latent heat`"                **surface-canopy_latent_heat**
   `"canopy sensible heat`"              **surface-canopy_sensible_heat**
   `"canopy surface water`"              **surface-canopy_surface_water**
   `"evapotranspiration`"                **surface-evapotranspiration**
   `"evaporation ground`"                **surface-evaporation_ground**
   `"evaporation litter`"                **surface-evaporation_litter**
   `"evaporation snow`"                  **surface-evaporation_snow**
   `"sublimation snow`"                  **surface-sublimation_snow**
   `"LAI`"                               **surface-LAI**
   `"SAI`"                               **surface-SAI**
   `"vegetation type`"                   **surface-vegetation_type**

   //General Flow Transport Energy
   `"mole fraction`"                     **mole_fraction**
   `"porosity`"                          **porosity**
   `"saturation liquid`"                 **saturation_liquid**
   `"saturation gas`"                    **saturation_gas**
   `"saturation ice`"                    **saturation_ice**
   `"water content`"                     **water_content**
   `"mass density liquid`"               **mass_density_liquid**
   `"mass density ice`"                  **mass_density_ice**
   `"mass density gas`"                  **mass_density_gas**
   `"density rock`"                      **density_rock**
   `"temperature`"                       **temperature**
   `"thermal conductivity`"              **thermal_conductivity**
   `"cell volume`"                       **cell_volume**

 */


#ifndef PKS_ECOSIM_HH_
#define PKS_ECOSIM_HH_

#include <map>
#include <vector>
#include <string>

#include "Epetra_MultiVector.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"

#include "VerboseObject.hh"
#include "TreeVector.hh"

#include "Key.hh"
#include "Mesh.hh"
#include "State.hh"
#include "BGCEngine.hh"
#include "PK_Factory.hh"
#include "PK_Physical_Default.hh"
#include "PK_Physical.hh"
#include "MeshPartition.hh"

namespace Amanzi {
namespace EcoSIM {

//using namespace Amanzi::Flow;

class EcoSIM : public PK_Physical_Default {

 public:

  //Unclear if the constructor is neccessary
  EcoSIM(Teuchos::ParameterList& pk_tree,
              const Teuchos::RCP<Teuchos::ParameterList>& plist,
              const Teuchos::RCP<State>& S,
              const Teuchos::RCP<TreeVector>& solution);
  // Virtual destructor
  ~EcoSIM();

  // is a PK
  // -- Setup data
  //virtual void Setup(const Teuchos::Ptr<State>&S);
  virtual void Setup() final;

  // -- initalize owned (dependent) variables
  //virtual void Initialize(const Teuchos::Ptr<State>& S);
  virtual void Initialize() final;

  // --provide timestep size
  virtual double get_dt() final {
    return dt_;
  }

  virtual void set_dt(double dt) final {
    dt_ = dt;
  }

  // -- commit the model
  //virtual void CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S);
  virtual void CommitStep(double t_old, double t_new, const Tag& tag) final;

  // -- Update diagnostics for vis.
  //virtual void CalculateDiagnostics(const Teuchos::RCP<State>& S) {}
  //virtual void CalculateDiagnostics(const Tag& tag) override;

  // -- advance the model
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit) final;

  virtual std::string name(){return "EcoSIM for ATS";};

  //This is not in the Alquimia_PK, for whatever reason it is defined in
  //The Chemistry_PK even though it isn't used there, and then included
  Teuchos::RCP<BGCEngine> bgc_engine() { return bgc_engine_; }

 private:

   //Helper functions from Alquimia
   void CopyToEcoSIM(int column,
           BGCProperties& props,
           BGCState& state,
           BGCAuxiliaryData& aux_data,
         const Tag& water_tag = Tags::DEFAULT);

   void CopyFromEcoSIM(const int cell,
                const BGCProperties& props,
                const BGCState& state,
                const BGCAuxiliaryData& aux_data,
              const Tag& water_tag = Tags::DEFAULT);

    //Helper functions from Alquimia
    void CopyToEcoSIM_process(int proc,
            BGCProperties& props,
            BGCState& state,
            BGCAuxiliaryData& aux_data,
          const Tag& water_tag = Tags::DEFAULT);

    void CopyFromEcoSIM_process(const int proc,
                 const BGCProperties& props,
                 const BGCState& state,
                 const BGCAuxiliaryData& aux_data,
               const Tag& water_tag = Tags::DEFAULT);

   int InitializeSingleProcess(int proc);

   int AdvanceSingleProcess(double dt, int proc);

   void ComputeNextTimeStep();

 protected:
  double dt_;
  double c_m_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_surf_; //might need this?
  Key domain_surface_;
  std::string passwd_ = "state";

  //The helper functions from BGC are protected not private (unclear why)
  //I don't think I need this here, probably in the engine
  void FieldToColumn_(AmanziMesh::Entity_ID column, const Epetra_Vector& vec,
          Teuchos::Ptr<Epetra_SerialDenseVector> col_vec);

  void FieldToColumn_(AmanziMesh::Entity_ID column, Teuchos::Ptr<Epetra_SerialDenseVector> vec,
          Teuchos::Ptr<Epetra_SerialDenseVector> col_vec);

  void MatrixFieldToColumn_(AmanziMesh::Entity_ID column, const Epetra_MultiVector& m_arr,
      Teuchos::Ptr<Epetra_SerialDenseMatrix> col_arr);

  //void FieldToColumn_(AmanziMesh::Entity_ID column, const Epetra_Vector& vec,
  //                                       double* col_vec);
  void ColDepthDz_(AmanziMesh::Entity_ID column,
                              Teuchos::Ptr<Epetra_SerialDenseVector> depth,
                              Teuchos::Ptr<Epetra_SerialDenseVector> dz);

  void VolDepthDz_(AmanziMesh::Entity_ID column,
                              Teuchos::Ptr<Epetra_SerialDenseVector> depth,
                              Teuchos::Ptr<Epetra_SerialDenseVector> dz,
			      Teuchos::Ptr<Epetra_SerialDenseVector> volume);

  void ColumnToField_(AmanziMesh::Entity_ID column, Epetra_Vector& vec,
                                 Teuchos::Ptr<Epetra_SerialDenseVector> col_vec);

  void ColumnToField_(AmanziMesh::Entity_ID column, Teuchos::Ptr<Epetra_SerialDenseVector> vec,
                                 Teuchos::Ptr<Epetra_SerialDenseVector> col_vec);

  void MatrixColumnToField_(AmanziMesh::Entity_ID column, Epetra_MultiVector& m_arr,
                                 Teuchos::Ptr<Epetra_SerialDenseMatrix> col_arr);

  int number_aqueous_components_;
  int ncells_per_col_;
  int num_columns_;
  int num_columns_local;
  int num_columns_global;
  int num_columns_global_ptype;
  int day0_, year0_, curr_day_, curr_year_;
  double saved_time_;
  double current_time_;
  double t_ecosim = 0.0;

  // keys
  Key mole_fraction_key_;
  Key porosity_key_;
  Key saturation_liquid_key_;
  Key saturation_gas_key_;
  Key saturation_ice_key_;
  Key elev_key_;
  Key water_content_key_;
  Key relative_permeability_key_;
  Key liquid_density_key_;
  Key ice_density_key_;
  Key gas_density_key_;
  Key gas_density_key_test_;
  Key rock_density_key_;
  Key T_key_;
  Key thermal_conductivity_key_;
  Key cell_volume_key_;
  Key ecosim_aux_data_key_;
  Key bulk_density_key_;
  Key hydraulic_conductivity_key_;
  Key sw_key_;
  Key lw_key_;
  Key air_temp_key_;
  Key vp_air_key_;
  Key wind_speed_key_;
  Key p_rain_key_;
  Key p_snow_key_;
  Key p_total_key_;
  Key f_wp_key_;
  Key f_root_key_;
  Key matric_pressure_key_;
  Key aspect_key_;
  Key slope_key_;
  Key lai_key_;
  Key sai_key_;
  Key v_type_key_;
  Key surface_energy_source_key_;
  Key subsurface_energy_source_key_;
  Key surface_water_source_key_;
  Key subsurface_water_source_key_;
  Key surface_energy_source_ecosim_key_;
  Key surface_water_source_ecosim_key_;
  Key subsurface_energy_source_ecosim_key_;
  Key subsurface_water_source_ecosim_key_;
  Key snow_depth_key_;
  Key snow_albedo_key_;
  Key canopy_lw_key_;
  Key canopy_latent_heat_key_;
  Key canopy_sensible_heat_key_;
  Key canopy_surface_water_key_;
  Key transpiration_key_;
  Key evaporation_canopy_key_;
  Key evaporation_ground_key_;
  Key evaporation_litter_key_;
  Key evaporation_snow_key_;
  Key sublimation_snow_key_;
  Key snow_temperature_key_;
  Key cap_pres_key_;

  Teuchos::RCP<BGCEngine> bgc_engine_;

  double atm_n2_, atm_o2_, atm_co2_, atm_ch4_, atm_n2o_, atm_h2_, atm_nh3_;
  double pressure_at_field_capacity, pressure_at_wilting_point;

 private:
  BGCState bgc_state_;
  BGCProperties bgc_props_;
  BGCAuxiliaryData bgc_aux_data_;
  BGCSizes bgc_sizes_;

  Teuchos::RCP<Epetra_SerialDenseVector> column_vol_save;
  Teuchos::RCP<Epetra_SerialDenseVector> column_wc_save;

  bool bgc_initialized_;
  bool has_energy, has_gas, has_ice, p_bool, a_bool, pheno_bool;
  std::vector<std::string> component_names_;
  int num_components;

 private:
  //factory registration
  static RegisteredPKFactory<EcoSIM> reg_;
};

} // namespace EcoSIM
} // namespace Amanzi

#endif
