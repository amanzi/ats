/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ATS version) (ecoon@lanl.gov)
*/

/*
  Delegate for modifying the predictor in the case of infiltration into dry soil.

  NOTE this uses only a domain, and assumes standard variable names.

*/

#ifndef PREDICTOR_DELEGATE_BC_FLUX_
#define PREDICTOR_DELEGATE_BC_FLUX_

#include "Mesh.hh"
#include "State.hh"

#include "TreeVector.hh"
#include "PDE_Diffusion.hh"
#include "wrm_partition.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {

class PredictorDelegateBCFlux {
 public:
  PredictorDelegateBCFlux(const Teuchos::RCP<const State>& S_next,
                          const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                          const Teuchos::RCP<Operators::PDE_Diffusion>& matrix,
                          const Teuchos::RCP<Flow::WRMPartition>& wrms,
                          std::vector<int>* bc_markers,
                          std::vector<double>* bc_values)
    : S_next_(S_next),
      mesh_(mesh),
      matrix_(matrix),
      wrms_(wrms),
      bc_markers_(bc_markers),
      bc_values_(bc_values)
  {}

  bool ModifyPredictor(double h, Teuchos::RCP<TreeVector> u)
  {
    return ModifyPredictor(u->Data().ptr());
  }

  bool ModifyPredictor(const Teuchos::Ptr<CompositeVector>& u);

 protected:
  class FluxBCFunctor {
   public:
    FluxBCFunctor(const Teuchos::RCP<std::vector<double>> Aff,
                  const Teuchos::RCP<std::vector<double>> lambda,
                  int face_index,
                  double cell_p,
                  double bc_flux,
                  double g_flux,
                  int dir,
                  double patm,
                  const Teuchos::RCP<Flow::WRM>& wrm)
      : Aff_(Aff),
        lambda_(lambda),
        face_index_(face_index),
        cell_p_(cell_p),
        bc_flux_(bc_flux),
        g_flux_(g_flux),
        wrm_(wrm),
        dir_(dir),
        patm_(patm)
    {}

    double operator()(double face_p) const
    {
      (*lambda_)[face_index_] = face_p;
      double s = wrm_->saturation(patm_ - face_p);
      double Krel = wrm_->k_relative(s);

      double q = 0;
      for (unsigned int n = 0; n != lambda_->size(); ++n)
        q += (*Aff_)[n] * (cell_p_ - (*lambda_)[n]);

      return (q + g_flux_) * Krel - bc_flux_;
    }

   protected:
   protected:
    Teuchos::RCP<std::vector<double>> Aff_;
    Teuchos::RCP<std::vector<double>> lambda_;
    int face_index_;
    double face_Mff_;
    double cell_p_;
    double bc_flux_;
    double g_flux_;
    int dir_;
    double patm_;
    Teuchos::RCP<Flow::WRM> wrm_;
  };

  struct Tol_ {
    Tol_(double eps)
      : eps_(eps)
    {}
    bool operator()(const double& a, const double& b) const { return std::abs(a - b) <= eps_; }
    double eps_;
  };

 protected:
  Teuchos::RCP<FluxBCFunctor> CreateFunctor_(int f,
                                             const Teuchos::Ptr<const CompositeVector>& pres);
  int CalculateLambda_(int f, const Teuchos::Ptr<const CompositeVector>& pres, double& lambda);

 protected:
  Teuchos::RCP<const State> S_next_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<Operators::PDE_Diffusion> matrix_;
  Teuchos::RCP<Flow::WRMPartition> wrms_;

  std::vector<int>* bc_markers_;
  std::vector<double>* bc_values_;
};

} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi

#endif
