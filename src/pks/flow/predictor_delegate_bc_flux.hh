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

  NOTE: this only works on HOST for now!
*/

#ifndef PREDICTOR_DELEGATE_BC_FLUX_
#define PREDICTOR_DELEGATE_BC_FLUX_

#include "Mesh.hh"
#include "State.hh"

#include "TreeVector.hh"
#include "PDE_Diffusion.hh"

#include "EvaluatorModelCVByMaterial.hh"
#include "wrm_model.hh"
#include "wrm_van_genuchten.hh"

namespace Amanzi {
namespace Flow {

class PredictorDelegateBCFlux {
 public:
  using WRMEval_type = EvaluatorModelCVByMaterial<Relations::WRMVanGenuchtenModel>;
  using WRMModel_type = WRMEval_type::Model_type;

  PredictorDelegateBCFlux(
    const Teuchos::RCP<const State>& S_next,
    const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
    const Teuchos::RCP<Operators::PDE_Diffusion>& matrix,
    const std::vector<std::pair<std::string, Teuchos::RCP<WRMModel_type>>>& wrms,
    const Teuchos::RCP<const CompositeVector_<int>> bc_markers,
    const Teuchos::RCP<const CompositeVector> bc_values)
    : S_next_(S_next),
      mesh_(mesh),
      matrix_(matrix),
      wrms_(wrms),
      bc_markers_(bc_markers),
      bc_values_(bc_values)
  {}

  bool ModifyPredictor(double h, Teuchos::RCP<TreeVector> u)
  {
    return ModifyPredictor(u->getData().ptr());
  }

  bool ModifyPredictor(const Teuchos::Ptr<CompositeVector>& u);

 protected:
  class FluxBCFunctor {
   public:
    FluxBCFunctor(const Kokkos::View<double*, Amanzi::DefaultHost>& Aff,
                  const Kokkos::View<double*, Amanzi::DefaultHost>& lambda,
                  int face_index,
                  double cell_p,
                  double bc_flux,
                  double g_flux,
                  int dir,
                  double patm,
                  const Relations::WRMVanGenuchten& wrm)
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

    // NOTE: this is not actually valid on device, but is decorated anyway to
    // silence warnings.
    KOKKOS_INLINE_FUNCTION
    double operator()(double face_p) const
    {
      lambda_(face_index_) = face_p;
      double s = wrm_.saturation(patm_ - face_p);
      double Krel = wrm_.k_relative(s);

      double q = 0;
      for (unsigned int n = 0; n != lambda_.extent(0); ++n) q += Aff_(n) * (cell_p_ - lambda_(n));

      return (q + g_flux_) * Krel - bc_flux_;
    }

   protected:
   protected:
    const Kokkos::View<double*, Amanzi::DefaultHost> Aff_, lambda_;
    int face_index_;
    double face_Mff_;
    double cell_p_;
    double bc_flux_;
    double g_flux_;
    int dir_;
    double patm_;
    const Relations::WRMVanGenuchten& wrm_;
  };

 protected:
  Teuchos::RCP<FluxBCFunctor> CreateFunctor_(AmanziMesh::Entity_ID f,
                                             AmanziMesh::Entity_ID c,
                                             const Teuchos::RCP<WRMModel_type>& wrm,
                                             const Teuchos::Ptr<const CompositeVector>& pres);

  int CalculateLambda_(AmanziMesh::Entity_ID f,
                       AmanziMesh::Entity_ID c,
                       const Teuchos::RCP<WRMModel_type>& wrm,
                       const Teuchos::Ptr<const CompositeVector>& pres,
                       double& lambda);

 protected:
  Teuchos::RCP<const State> S_next_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<Operators::PDE_Diffusion> matrix_;
  std::vector<std::pair<std::string, Teuchos::RCP<WRMModel_type>>> wrms_;
  Teuchos::RCP<const CompositeVector_<int>> bc_markers_;
  Teuchos::RCP<const CompositeVector> bc_values_;
};

} // namespace Flow
} // namespace Amanzi

#endif
