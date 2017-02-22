/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
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
#include "Epetra_SerialDenseVector.h"

#include "VerboseObject.hh"
#include "TreeVector.hh"

#include "pk_factory.hh"
#include "pk_physical_base.hh"

#include "SoilCarbonParameters.hh"
#include "PFT.hh"
#include "SoilCarbon.hh"

namespace Amanzi {
namespace BGC {

class BGCSimple : public PKPhysicalBase {

 public:
  BGCSimple(const Teuchos::RCP<Teuchos::ParameterList>& plist,
            Teuchos::ParameterList& FElist,
            const Teuchos::RCP<TreeVector>& solution);

  // is a PK
  // -- Setup data
  virtual void setup(const Teuchos::Ptr<State>& S);

  // -- Initialize owned (dependent) variables.
  virtual void initialize(const Teuchos::Ptr<State>& S);

  // -- provide a timestep size
  virtual double get_dt() {
    return dt_;
  }

  // -- Commit any secondary (dependent) variables.
  virtual void commit_state(double dt, const Teuchos::RCP<State>& S);

  // -- Update diagnostics for vis.
  virtual void calculate_diagnostics(const Teuchos::RCP<State>& S) {}

  // -- advance the model
  virtual bool advance(double dt);

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
      mesh.face_get_cells(col_face, AmanziMesh::USED, &facecells);
      ASSERT(facecells.size() == 1);

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
  Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh_;

  // physical structs needed by model
  std::vector<Teuchos::RCP<SoilCarbonParameters> > sc_params_;
  std::vector<std::vector<Teuchos::RCP<PFT> > > pfts_;       // this also contains state data!
  std::vector<std::vector<Teuchos::RCP<PFT> > > pfts_old_;   // need two copies for failed timesteps
  std::vector<std::vector<Teuchos::RCP<SoilCarbon> > > soil_carbon_pools_;

  // evaluator for transpiration
  Teuchos::RCP<PrimaryVariableFieldEvaluator> trans_eval_;
  Teuchos::RCP<PrimaryVariableFieldEvaluator> sw_eval_;
  Teuchos::RCP<PrimaryVariableFieldEvaluator> lai_eval_;
  
  // keys
  Key shortwave_incoming_key_;
  Key shortwave_shaded_key_;

  // extras
  double lat_;
  double wind_speed_ref_ht_;
  double cryoturbation_coef_;
  int ncells_per_col_;
  std::string soil_part_name_;

 private:
  // factory registration
  static RegisteredPKFactory<BGCSimple> reg_;


};

} // namespace
} // namespace


#endif
