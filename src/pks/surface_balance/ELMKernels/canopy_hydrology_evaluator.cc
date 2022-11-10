
#include "Key.hh"

namespace Amanzi {
namespace ELMKernels {

CanopyHydrologyEvaluator::CanopyHydrologyEvaluator(Teuchos::ParameterList& plist) :
    EvaluatorSecondaryMonotypeCV(plist)
{
  // determine the domain
  Key akey = my_keys_.front().first;
  domain_ = Keys::getDomain(akey);
  akey = Keys::getVarName(akey);
  domain_snow_ = Keys::readDomainHint(plist_, domain_, "surface", "snow");
  auto tag = my_keys_.front().second;

  // my keys
  // -- sources
  // gather all dependencies, assign any default parameters
  // would prefer to do this as simply as possible, 
  // because the function signatures will change soon
}

// Required methods from EvaluatorSecondaryMonotypeCV
void
CanopyHydrologyEvaluator::Evaluate_(const State& S,
        const std::vector<CompositeVector*>& results)
{
  auto mesh = S.GetMesh(domain_);
  auto tag = my_keys_.front().second;

  // collect dependencies

  // collect output vecs

  // call canopy_hydrology kernels
  ELM::kokkos_canopy_hydrology(*S, atm_forcing->forc_PREC, dtime, time_plus_half_dt)

}

void
CanopyHydrologyEvaluator::EvaluatePartialDerivative_(const State& S,
        const Key& wrt_key, const Tag& wrt_tag,
        const std::vector<CompositeVector*>& results) {//nope}


// custom EC used to set subfield names
void
CanopyHydrologyEvaluator::EnsureCompatibility_Structure_(State& S)
{
  S.GetRecordSetW(my_keys_.front().first).set_subfieldnames(
    {"bare", "water", "snow"});
}


void CanopyHydrologyEvaluator::EnsureCompatibility_ToDeps_(State& S)
{
  // new state!
  if (land_cover_.size() == 0)
    land_cover_ = getLandCover(S.ICList().sublist("land cover types"),
            {"albedo_ground", "emissivity_ground"});

  for (auto dep : dependencies_) {
    auto& fac = S.Require<CompositeVector,CompositeVectorSpace>(dep.first, dep.second);
    if (Keys::getDomain(dep.first) == domain_snow_) {
      fac.SetMesh(S.GetMesh(domain_snow_))
        ->SetGhosted()
        ->AddComponent("cell", AmanziMesh::CELL, 1);
    } else {
      fac.SetMesh(S.GetMesh(domain_))
        ->SetGhosted()
        ->AddComponent("cell", AmanziMesh::CELL, 1);
    }
  }
}

}  // namespace Relations
}  // namespace ELMKernels
