/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! A simple conservation equation.

/*!

This is a helper PK that does not provide sufficient information to be a
standalone PK, but provides tools for conservation equations.  This includes
dealing with all needed keys, data, and evaluators for the conserved quantity
and source terms.

.. _conservation-equation-pk-spec:
.. admonition:: conservation-equation-pk-spec

    * `"domain`" ``[string]`` Mesh on which the balance is to be done.

    * `"conserved quantity key`" ``[string]`` The conserved quantity.

    * `"is source term`" ``[bool]`` **true**

    * `"source key`" ``[string]`` **DOMAIN-source_sink** Units are in conserved
      quantity per second per cell volume, :math:`Q`.

    * `"source term is differentiable`" ``[bool]`` **true**

    * `"finite difference source term`" ``[bool]`` **false** If the source is
      not differentiable, should we attempt to form a local finite difference
      to calculate dQ/du?
    
    * `"time discretization theta`" ``[double]`` **1.0** :math:`\theta` in a
      Crank-Nicholson time integration scheme.  1.0 implies fully implicit, 0.0
      implies explicit, 0.5 implies C-N.  Note, only used in the implicit
      scheme, and may only be used for the source term.

    * `"modify predictor positivity preserving`" ``[bool]`` **false** If true,
      predictors are modified to ensure that the primary variable is always >
      0.  These systems may be stiff, and this does not guarantee positivity,
      so time integration methods may need to be chosen with care.

    * `"absolute error tolerance`" ``[double]`` **550.0** a_tol in the standard
      error norm calculation.  Defaults to a small amount of water.  Units are
      the same as the conserved quantity.

*/

#pragma once

#include <cmath>


#include "Key.hh"
#include "Reductions.hh"
#include "TreeVector.hh"
#include "SolverDefs.hh"

#include "PK_Adaptors.hh"
#include "PK_MixinImplicit.hh"
#include "PK_MixinExplicit.hh"
#include "PK_MixinLeaf.hh"
#include "PK_Default.hh"
#include "PK_Factory.hh"


namespace ATS {
namespace Basic {

using namespace Amanzi;

template <class Base_t>
class PK_MixinConservationEquation : public Base_t {

public:

  PK_MixinConservationEquation(const Teuchos::RCP<Teuchos::ParameterList>& pk_tree,
              const Teuchos::RCP<Teuchos::ParameterList>& global_plist,
              const Teuchos::RCP<State>& S)
    : Base_t(pk_tree, global_plist, S)
  {
    layer_ = plist_->template get<std::string>("layer name", domain_);

    // process keys
    conserved_quantity_key_ = Keys::readKey(*plist_, layer_, "conserved quantity");

    is_source_ = plist_->template get<bool>("is source term", true);
    if (is_source_) {
      source_is_extensive_ = plist_->template get<bool>("source is extensive", false);
      source_key_ = Keys::readKey(*plist_, layer_, "source", Keys::getVarName(conserved_quantity_key_)+"_source");
      if (!source_is_extensive_) {
        intensive_source_key_ = source_key_;
        source_key_ = intensive_source_key_ + "_times_volume";
      }
    }

    cv_key_ = Keys::readKey(*plist_, domain_, "cell volume", "cell_volume");
    
    // process error norm constants
    atol_ = plist_->template get<double>("absolute error tolerance factor", 1.);
    rtol_ = plist_->template get<double>("relative error tolerance factor", 1.);
    fluxtol_ = plist_->template get<double>("flux error tolerance factor", 1.);

    // process globalization
    modify_predictor_positivity_preserving_ = plist_->template get<bool>("modify predictor positivity preserving", false);
    modify_correction_positivity_preserving_ = plist_->template get<bool>("modify correction positivity preserving", false);
    admit_only_positivity_preserving_ = plist_->template get<bool>("admit only positivity preserving", false);
  }

  // -- Setup data.
  void SetupAtTag(const Key& tag) {
    // require primary evaluator and set the structure of that field
    S_->template Require<CompositeVector, CompositeVectorSpace>(key_, tag, key_)
        .SetMesh(mesh_)
        ->SetComponent("cell", AmanziMesh::CELL, 1);
    S_->RequireEvaluator(key_, tag);

    // similarly set the structure of the conserved quantity
    S_->template Require<CompositeVector, CompositeVectorSpace>(conserved_quantity_key_, tag)
        .SetMesh(mesh_)
        ->AddComponent("cell", AmanziMesh::CELL, 1);

    if (is_source_) {
      // set the structure of the source
      S_->template Require<CompositeVector, CompositeVectorSpace>(source_key_, tag)
          .SetMesh(mesh_)
          ->AddComponent("cell", AmanziMesh::CELL, 1);

      // add an evaluator for converting the default, an intensive quantity, into an extensive one.
      if (!source_is_extensive_) {
        // source times volume
        Teuchos::ParameterList& source_times_vol = S_->FEList().sublist(source_key_);
        if (!source_times_vol.isParameter("tag"))
          source_times_vol.set("tag", tag);
        if (!source_times_vol.isParameter("evaluator type"))
          source_times_vol.set("evaluator type", "multiplicative");
        if (!source_times_vol.isParameter("dependencies")) {
          std::vector<std::string> deps{intensive_source_key_, cv_key_};
          source_times_vol.set("dependencies", Teuchos::Array<std::string>(deps));
        }
      }
      S_->RequireEvaluator(source_key_, tag);
    }
  }

  double
  ErrorNorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du) {
    using Reduction_t = Amanzi::Reductions::MaxLocArray<int,double,GO>;
    using Reductor_t = Amanzi::Reductions::MaxLoc<double, GO>;
    
    
    // Abs tol based on old conserved quantity -- we know these have been vetted
    // at some level whereas the new quantity is some iterate, and may be
    // anything from negative to overflow.
    S_->GetEvaluator(conserved_quantity_key_).Update(*S_, this->name());
    const auto& conserved = S_->template Get<CompositeVector>(conserved_quantity_key_, tag_old_);
    const auto& cv = S_->template Get<CompositeVector>(cv_key_, tag_old_);

    // VerboseObject stuff.
    Teuchos::OSTab tab = vo_->getOSTab();
    if (vo_->os_OK(Teuchos::VERB_MEDIUM))
      *vo_->os() << "ENorm (Infnorm) of: " << conserved_quantity_key_ << ": " << std::endl;

    double dt = S_->time(tag_new_) - S_->time(tag_old_);
    Reductor_t enorm_all;

    for (const auto& comp : *du->Data()) {
      Reductor_t enorm_comp;

      auto du_v = du->Data()->template ViewComponent<>(comp,false);
      AMANZI_ASSERT(du_v.extent(1) == 1);
      
      if (comp == "cell") {
        // error done relative to extensive, conserved quantity
        auto cv_v = cv.template ViewComponent<>(comp, false);
        auto conserved_v = conserved.template ViewComponent<>(comp, false);
        
        Kokkos::parallel_reduce("ConservationEquation::ErrorNorm(cell) reduction",
                du_v.extent(0),
                KOKKOS_LAMBDA(const int& c, Reductor_t& enorm) {
                  double local_enorm = abs(dt * du_v(c,0)) / cv_v(c,0) / (atol_ + rtol_ * abs(conserved_v(c,0)));
                  //if (c == 99) std::cout << std::setprecision(16) << "local_enorm = " << dt << " * " << du_v(c,0) << " / " << cv_v(c,0) << " / " << atol_ << " + " << rtol_ << " * " << conserved_v(c,0) << " = " << local_enorm << std::endl;
                  Reductor_t l_enorm(local_enorm, c);
                  enorm += l_enorm;
                }, enorm_comp);

      } else if (comp == "face") {
        // error in flux -- relative to cell's extensive conserved quantity
        cv.ScatterMasterToGhosted("cell");
        auto cv_v = cv.template ViewComponent<>("cell", false);
        conserved.ScatterMasterToGhosted();
        auto conserved_v = conserved.template ViewComponent<>("cell", false);

        const AmanziMesh::Mesh* m = mesh_.get();
        
        Kokkos::parallel_reduce("ConservationEquation::ErrorNorm(face) reduction",
                du_v.extent(0),
                KOKKOS_LAMBDA(const int& f, Reductor_t& enorm) {
                  AmanziMesh::Entity_ID_View cells;
                  m->face_get_cells(f, AmanziMesh::Parallel_type::OWNED, cells);
                  double cv_min = cells.extent(0) == 1 ? cv_v(cells(0),0) :
                                  fmin(cv_v(cells(0),0), cv_v(cells(1),1));
                  double conserved_min = cells.extent(0) == 1 ? conserved_v(cells(0),0) :
                                  fmin(conserved_v(cells(0),0), conserved_v(cells(1),1));

                  double local_enorm = fluxtol_ * dt * abs(du_v(f,0)) / cv_min / (atol_ + rtol_ * abs(conserved_min));
                  Reductor_t l_enorm(local_enorm, f);
                  enorm += l_enorm;
                }, enorm_comp);
      }
      
      // note, now it is GID
      enorm_comp.val.second = du->Data()->getMap()->ComponentMap(comp)
                       ->getGlobalElement(enorm_comp.val.second);

      // Write out Inf norms too.
      if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
        Teuchos::Array<double> infnorm(1);
        du->Data()->GetComponent(comp, false)->normInf(infnorm);

        Reductor_t err;
        Teuchos::reduceAll(*du->Data()->Comm(), Reduction_t(), 1, &enorm_comp.val, &err.val);
        *vo_->os() << "  ENorm (" << comp << ") = " << err.val.first << "[" << err.val.second << "] (" << infnorm[0] << ")" << std::endl;
      }

      enorm_all += enorm_comp;
    }

    Reductor_t err;
    Teuchos::reduceAll(*du->Data()->Comm(), Reduction_t(), 1, &enorm_all.val, &err.val);
    return enorm_all.val.first;
  }

  void ChangedSolution() {
    this->ChangedSolutionPK(tag_new_);
  }

  // -- limit changes in a valid time step
  bool ValidStep(const Key& tag_old, const Key& tag_new) {
    bool valid =  Base_t::ValidStep(tag_old, tag_new);
    if (!valid) return valid;
    
    if (admit_only_positivity_preserving_) {
      Errors::Message msg("ConservationEquation: \"admit only positivity preserving\" is not yet implemented.");
      throw(msg);
    }
    return valid;
  }

  bool ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0,
                       Teuchos::RCP<TreeVector> u) {
    bool result = Base_t::ModifyPredictor(h, u0, u);
    if (modify_predictor_positivity_preserving_) {
      Errors::Message msg("ConservationEquation: \"modify predictor positivity preserving\" is not yet implemented.");
      throw(msg);
    }
    return result;
  }

  // -- Possibly modify the correction before it is applied
  AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
  ModifyCorrection(double h, Teuchos::RCP<const TreeVector> res,
                   Teuchos::RCP<const TreeVector> u,
                   Teuchos::RCP<TreeVector> du){
    auto base_changed = Base_t::ModifyCorrection(h, res, u, du);
    if (modify_correction_positivity_preserving_) {
      auto u_v = u->Data()->ViewComponent("cell", false);
      auto du_v = du->Data()->ViewComponent("cell", false);
      int changed = 0;
      Kokkos::parallel_reduce(
          "PK_MixinConservationEquation::ModifyCorrection preserving positivity",
          du_v.extent(0),
          KOKKOS_LAMBDA(const int& i, int& count) {
            if (u_v(i,0) - du_v(i,0) < 0.) {
              std::cout << "POS PRESERV: " << i << "u = " << u_v(i,0) << " du = " << du_v(i,0) << std::endl;
              du_v(i,0) = u_v(i,0);
              count += 1;
            }
          }, changed);
      if (changed) return AmanziSolvers::FnBaseDefs::CORRECTION_MODIFIED_LAG_BACKTRACKING;
    }
    return base_changed;
  }

 protected:
  using Base_t::tag_new_;
  using Base_t::tag_old_;
  using Base_t::S_;
  using Base_t::key_;
  using Base_t::mesh_;
  using Base_t::db_;
  using Base_t::plist_;
  using Base_t::domain_;
  using Base_t::vo_;

  Key layer_;
  Key conserved_quantity_key_;
  Key cv_key_;
  Key source_key_;
  Key intensive_source_key_;
  
  bool is_source_, source_is_extensive_;

  double theta_;

  bool modify_predictor_positivity_preserving_;
  bool modify_correction_positivity_preserving_;
  bool admit_only_positivity_preserving_;

  double atol_, rtol_, fluxtol_;
  
};

}  // namespace Basic
}  // namespace Amanzi



