/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/
//! Calculates a face value from cell values through upwinding.

#include "EvaluatorUpwindPotentialDifference.hh"

namespace Amanzi {

EvaluatorUpwindPotentialDifference::EvaluatorUpwindPotentialDifference(
    Teuchos::ParameterList& plist)
    : EvaluatorSecondary(plist) {
  if (dependencies_.size() == 0) {
    // get the cell-based coefficient
    Key my_tag = my_keys_[0].second;
    Key my_key = my_keys_[0].first;
    KeyPair domain_key = Keys::splitKey(my_key);

    if (Keys::startsWith(domain_key.second, "upwind_")) {
      Key my_default = domain_key.second.substr(7,domain_key.second.size());
      dependencies_.emplace_back(std::make_pair(
          Keys::readKey(plist, domain_key.first, "upwinded field", my_default),
          my_tag));
    } else {
      dependencies_.emplace_back(std::make_pair(
          Keys::readKey(plist, domain_key.first, "upwinded field"),
          my_tag));
    }

    // get the potential field
    Key potential_key = Keys::readKey(plist, domain_key.first, "potential field");
    dependencies_.emplace_back(std::make_pair(potential_key, my_tag));    
  }
}

EvaluatorUpwindPotentialDifference&
EvaluatorUpwindPotentialDifference::operator=(const Evaluator& other)
{
  auto other_as_eval = dynamic_cast<const EvaluatorUpwindPotentialDifference*>(&other);
  AMANZI_ASSERT(other_as_eval != nullptr);
  return this->operator=(*other_as_eval);
}

void
EvaluatorUpwindPotentialDifference::EnsureCompatibility(State& S)
{
  AMANZI_ASSERT(my_keys_.size() == 1);
  AMANZI_ASSERT(dependencies_.size() == 2);
  auto my_key = my_keys_[0];
  auto& my_fac = S.Require<CompositeVector,CompositeVectorSpace>(my_key.first, my_key.second, my_key.first);

  if (my_fac.Mesh().get()) {
    my_fac.SetComponent("face", AmanziMesh::FACE, 1)->SetGhosted();

    for (const auto& dep : dependencies_) {
      auto& dep_fac = S.Require<CompositeVector,CompositeVectorSpace>(dep.first, dep.second);
      dep_fac.SetMesh(my_fac.Mesh());
      dep_fac.AddComponent("cell", AmanziMesh::CELL, 1)->SetGhosted();
      S.RequireEvaluator(dep.first, dep.second).EnsureCompatibility(S);

      if (S.HasDerivativeSet(my_key.first, my_key.second)) {
        for (const auto& deriv : S.GetDerivativeSet(my_key.first, my_key.second)) {
          auto wrt = Keys::splitKeyTag(deriv.first);
          if (S.GetEvaluator(dep.first, dep.second).IsDifferentiableWRT(S, wrt.first, wrt.second)) {
            S.RequireDerivative<CompositeVector,CompositeVectorSpace>(
                dep.first, dep.second, wrt.first, wrt.second).Update(my_fac);
          }
        }
      }
    }
  }
  EnsureCompatibility_Flags_(S);
}

void
EvaluatorUpwindPotentialDifference::Update_(State& S)
{
  const auto& cells = S.Get<CompositeVector>(dependencies_[0].first, dependencies_[0].second);
  const auto& potential = S.Get<CompositeVector>(dependencies_[1].first, dependencies_[1].second);
  auto& faces = S.GetW<CompositeVector>(my_keys_[0].first, my_keys_[0].second, my_keys_[0].first);
  Impl::upwindCellToFace(cells, potential, faces);
}



void
EvaluatorUpwindPotentialDifference::UpdateDerivative_(State& S,
        const Key& wrt_key, const Key& wrt_tag)
{
  AMANZI_ASSERT(false);
}


bool
EvaluatorUpwindPotentialDifference::UpdateDerivative(State& S,
        const Key& requestor, const Key& wrt_key, const Key& wrt_tag)
{
  Key wrt = Keys::getKeyTag(wrt_key, wrt_tag);
  auto& deriv_request_set = deriv_requests_[wrt];
  bool update = deriv_request_set.empty(); // not done once

  AMANZI_ASSERT(S.GetEvaluator(dependencies_[0].first, dependencies_[0].second)
                .IsDifferentiableWRT(S, wrt_key, wrt_tag));

  Teuchos::OSTab tab = vo_.getOSTab();
  if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
    *vo_.os() << "Variable d" << my_keys_[0].first << "_d" << wrt_key
              << " requested by " << requestor << std::endl;
  }

  Key my_request = Key{ "d" } + Keys::getKeyTag(my_keys_[0].first, my_keys_[0].second) +
                                  "_d" + Keys::getKeyTag(wrt_key, wrt_tag);

  // differentiate the cell component and upwind that
  update |= S.GetEvaluator(dependencies_[0].first, dependencies_[0].second)
             .UpdateDerivative(S, my_request, wrt_key, wrt_tag);

  update |= S.GetEvaluator(dependencies_[1].first, dependencies_[1].second)
             .Update(S, my_request);

  if (update) {
    if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
      *vo_.os() << "  ... updating." << std::endl;
    }
    const auto& potential = S.Get<CompositeVector>(dependencies_[1].first,
            dependencies_[1].second); 
    const auto& dcell_dp = S.GetDerivative<CompositeVector>(dependencies_[0].first,
            dependencies_[0].second, wrt_key, wrt_tag);
    auto& dface_dp = S.GetDerivativeW<CompositeVector>(my_keys_[0].first,
            my_keys_[0].second, wrt_key, wrt_tag, my_keys_[0].first);
    Impl::upwindCellToFace(dcell_dp, potential, dface_dp);

    deriv_request_set.clear();
    deriv_request_set.insert(requestor);
    return true;

  } else {
    // Otherwise, simply service the request
    if (deriv_request_set.find(requestor) == deriv_request_set.end()) {
      if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
        *vo_.os() << "  ... not updating but new to this request." << std::endl;
      }
      deriv_request_set.insert(requestor);
      return true;
    } else {
      if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
        *vo_.os() << "  ... has not changed." << std::endl;
      }
      return false;
    }
  }
}


namespace Impl {   
void upwindCellToFace(const CompositeVector& cells,
        const CompositeVector& potential, CompositeVector& faces)
{
  cells.ScatterMasterToGhosted("cell");
  potential.ScatterMasterToGhosted("cell");
  
  { // scope for views
    const AmanziMesh::Mesh* m = cells.getMap()->Mesh().get();
    auto cells_v = cells.ViewComponent("cell");
    auto potential_v = potential.ViewComponent("cell");
    auto faces_v = faces.ViewComponent("face", false);

    Kokkos::parallel_for(
        "EvaluatorUpwindPotentialDifference",
        faces_v.extent(0),
        KOKKOS_LAMBDA(const int& f) {
          AmanziMesh::Entity_ID_View cells;
          m->face_get_cells(f, AmanziMesh::Parallel_type::ALL, cells);
          if (cells.size() == 1) {
            faces_v(f,0) = cells_v(cells(0),0);
          } else {
            faces_v(f,0) = potential_v(cells(0),0) > potential_v(cells(1),0) ?
                           cells_v(cells(0),0) :
                           (potential_v(cells(1),0) > potential_v(cells(0),0) ?
                            cells_v(cells(1),0) :
                            (cells_v(cells(0),0) + cells_v(cells(1),0)) / 2.0);
          }
        });
  }
}
} // namespace Impl
  
} // namespace Amanzi
