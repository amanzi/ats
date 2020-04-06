/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
           Chonggang Xu (cxu@lanl.gov)
*/
//! Above and below-ground carbon cycle model.

/*!

This is a multi-leaf layer, big-leaf vegetation model coupled to a Century
model for belowground carbon decomposition.

It leverages a PFT-based structure which allows multiple height-sorted PFTs to
coexist on the same grid cells, with the shorter PFTs getting whatever light is
left in the understory.

The implementation is based on an old, standalone code by Chonggang Xu, and
adapted for ATS.  While this is not simple, it is called BGC simple as it is
about the least amount of complexity required to get a reasonable carbon cycle
into ATS.

Outputs of this include transpiration, a critical sink for hydrology, as it
solves photosynthesis based on water availability.

Note this is an "explicit update PK," or effectively a forward Euler timestep
that is not written in ODE form.

Note this works on both the surface (vegetation) and subsurface (decomposition)
meshes.  **It is required** that the subsurface mesh is a "columnar" mesh, and
that build_columns in the subsurface Mesh_ spec has been supplied.

.. _bgc-simple-spec:
.. admonition:: bgc-simple-spec

  * `"initial time step`" ``[double]`` **1.0** Initial time step size `[s]`

  * `"number of carbon pools`" ``[int]`` **7** Unclear whether this can actually change?

  * `"soil carbon parameters`" ``[soil-carbon-spec-list]`` List of soil carbon parameters by soil mesh partition region name.

  * `"pft parameters`" ``[pft-spec-list]`` List of PFT parameters by PFT name.

  * `"latitude [degrees]`" ``[double]`` **60** Latitude of the simulation in degrees.  Used in radiation balance.

  * `"wind speed reference height [m]`" ``[double]`` **2.0** Reference height of the wind speed dataset.

  * `"cryoturbation mixing coefficient [cm^2/yr]`" ``[double]`` **5.0** Controls diffusion of carbon into the subsurface via cryoturbation.

  * `"leaf biomass initial condition`" ``[initial-conditions-spec]`` Sets the leaf biomass IC.

  * `"domain name`" ``[string]`` **domain**

  * `"surface domain name`" ``[string]`` **surface**

  * `"transpiration key`" ``[string]`` **DOMAIN-transpiration** The distributed transpiration flux `[mol s^-1]`

  * `"shaded shortwave radiation key`" ``[string]``
    **SURFACE_DOMAIN-shaded_shortwave_radiation** Shortwave radiation that gets
    past the canopy and teo the bare ground for soil evaporation. `[W m^-2]`

  * `"total leaf area index key`" ``[string]`` **SURFACE_DOMAIN-total_leaf_area_index** Total LAI across all PFTs.

  EVALUATORS:

  - `"temperature`" The soil temperature `[K]`
  - `"pressure`" soil mafic potential `[Pa]`
  - `"surface-cell_volume`" `[m^2]`
  - `"surface-incoming shortwave radiation`" `[W m^-2]`
  - `"surface-air_temperature`" `[K]`
  - `"surface-relative_humidity`" `[-]`
  - `"surface-wind_speed`" `[m s^-1]`
  - `"surface-co2_concentration`" `[ppm]`
  
*/

#ifndef PKS_BGC_SIMPLE_HH_
#define PKS_BGC_SIMPLE_HH_

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Epetra_SerialDenseVector.h"

#include "VerboseObject.hh"
#include "TreeVector.hh"

#include "PK_Factory.hh"
#include "pk_physical_default.hh"

#include "SoilCarbonParameters.hh"
#include "PFT.hh"
#include "SoilCarbon.hh"

namespace Amanzi {
namespace BGC {

class BGCSimple : public PK_Physical_Default {

 public:

  BGCSimple(Teuchos::ParameterList& pk_tree,
            const Teuchos::RCP<Teuchos::ParameterList>& glist,
            const Teuchos::RCP<State>& S,

            const Teuchos::RCP<TreeVector>& solution);

  // is a PK
  // -- Setup data
  virtual void Setup(const Teuchos::Ptr<State>& S);

  // -- Initialize owned (dependent) variables.
  virtual void Initialize(const Teuchos::Ptr<State>& S);

  // -- provide a timestep size
  virtual double get_dt() {
    return dt_;
  }
  virtual void set_dt(double dt) {
    dt_ = dt;
  }

  // -- Commit any secondary (dependent) variables.
  virtual void CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S);

  // -- Update diagnostics for vis.
  virtual void CalculateDiagnostics(const Teuchos::RCP<State>& S) {}

  // -- advance the model
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit);

  virtual std::string name(){return "BGC simple";};

 protected:
  void FieldToColumn_(AmanziMesh::Entity_ID col, const Epetra_Vector& vec,
                     Teuchos::Ptr<Epetra_SerialDenseVector> col_vec,
                     bool copy=true);
  void ColDepthDz_(AmanziMesh::Entity_ID col,
                   Teuchos::Ptr<Epetra_SerialDenseVector> depth,
                   Teuchos::Ptr<Epetra_SerialDenseVector> dz);

  class ColIterator {
   public:
    ColIterator(const AmanziMesh::Mesh& mesh,
                AmanziMesh::Entity_ID col_face, int ncells=0) {
      if (ncells > 0) cells_.reserve(ncells);
      AmanziMesh::Entity_ID_List facecells;
      mesh.face_get_cells(col_face, AmanziMesh::Parallel_type::ALL, &facecells);
      AMANZI_ASSERT(facecells.size() == 1);

      AmanziMesh::Entity_ID c = facecells[0];
      while (c >= 0) {
        cells_.push_back(c);
        c = mesh.cell_get_cell_below(c);
      }
    }

    typedef std::vector<AmanziMesh::Entity_ID>::const_iterator const_iterator;

    const_iterator begin() { return cells_.begin(); }
    const_iterator end() { return cells_.end(); }
    std::size_t size() { return cells_.size(); }
    AmanziMesh::Entity_ID operator[](std::size_t i) { return cells_[i]; }

   private:
    std::vector<AmanziMesh::Entity_ID> cells_;
  };


 protected:
  double dt_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_surf_;
  Key domain_surf_;
  
  // physical structs needed by model
  std::vector<Teuchos::RCP<SoilCarbonParameters> > sc_params_;
  std::vector<std::vector<Teuchos::RCP<PFT> > > pfts_;       // this also contains state data!
  std::vector<std::vector<Teuchos::RCP<PFT> > > pfts_old_;   // need two copies for failed timesteps
  std::vector<std::vector<Teuchos::RCP<SoilCarbon> > > soil_carbon_pools_;

  // evaluator for transpiration
  Teuchos::RCP<PrimaryVariableFieldEvaluator> trans_eval_;
  Teuchos::RCP<PrimaryVariableFieldEvaluator> sw_eval_;
  Teuchos::RCP<PrimaryVariableFieldEvaluator> lai_eval_;
  
  // extras
  double lat_;
  double wind_speed_ref_ht_;
  double cryoturbation_coef_;
  int ncells_per_col_;
  std::string soil_part_name_;

  // keys
  Key trans_key_;
  Key shaded_sw_key_;
  Key total_lai_key_;
  
 private:
  // factory registration
  static RegisteredPKFactory<BGCSimple> reg_;


};

} // namespace
} // namespace


#endif
