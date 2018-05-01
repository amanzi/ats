/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  An evaluator for pulling the darcy flux, at the surface, from the
  subsurface field and putting it into a surface field.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_RELATIONS_OVERLAND_SOURCE_FROM_SUBSURFACE_FLUX_EVALUATOR_HH_
#define AMANZI_RELATIONS_OVERLAND_SOURCE_FROM_SUBSURFACE_FLUX_EVALUATOR_HH_

#include "Evaluator_Factory.hh"
#include "EvaluatorSecondary.hh"

namespace Amanzi {
namespace Relations {

class OverlandSourceFromSubsurfaceFluxEvaluator :
    public EvaluatorSecondary<CompositeVector, CompositeVectorSpace> {

 public:
  explicit
  OverlandSourceFromSubsurfaceFluxEvaluator(Teuchos::ParameterList& plist);

  OverlandSourceFromSubsurfaceFluxEvaluator(const OverlandSourceFromSubsurfaceFluxEvaluator& other);

  Teuchos::RCP<Evaluator> Clone() const;

  // custom ensure compatibility as all data is not just on the same components
  virtual void EnsureCompatibility(State& S);

protected:
  // Required methods from SecondaryVariableFieldEvaluator
  virtual void Evaluate_(const State& S,
                         CompositeVector& result);
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key, const Key& wrt_tag,
                                          CompositeVector& result);

  void IdentifyFaceAndDirection_(const State& S);

  typedef std::pair<int, double> FaceDir;
  Teuchos::RCP<std::vector<FaceDir> > face_and_dirs_;

  Key flux_key_;
  Key dens_key_;
  bool volume_basis_;

  Key surface_mesh_key_;
  Key subsurface_mesh_key_;

 private:
  static Utils::RegisteredFactory<Evaluator,OverlandSourceFromSubsurfaceFluxEvaluator> fac_;

};

} //namespace
} //namespace

#endif
