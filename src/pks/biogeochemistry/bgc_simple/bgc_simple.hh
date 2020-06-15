/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon, Chonggang Xu

   Simple implementation of CLM's Century model for carbon decomposition and a
   simplified 2-PFT (sedge, moss) vegetation model for creating carbon.

   ------------------------------------------------------------------------- */

#ifndef PKS_BGC_SIMPLE_HH_
#define PKS_BGC_SIMPLE_HH_

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

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
