/*
  Copyright 2010-201x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (coonet@ornl.gov)
*/

//! Water content from ponded depth.

/*!

*/

#pragma once

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "Key.hh"
#include "Factory.hh"
#include "StateDefs.hh"
#include "EvaluatorModel_CompositeVector.hh"

namespace ATS {
namespace Flow {
namespace Relations {

using namespace Amanzi;


template <class InView_type, class OutView_type>
class SurfaceWaterContent {
 public:
  static const int n_dependencies = 2;
  static const std::string name;

  SurfaceWaterContent(Teuchos::ParameterList& plist)
  {
    wc_key_ = Keys::cleanPListName(plist.name());
    std::string domain = Keys::getDomain(wc_key_);
    pd_key_ = Keys::readKey(plist, domain, "surface water depth", "ponded_depth");
    dens_key_ = Keys::readKey(plist, domain, "density", "molar_density_liquid");
  }

  void SetViews(const std::vector<InView_type>& dependency_views,
                const std::vector<OutView_type>& result_views,
                const State& S)
  {
    AMANZI_ASSERT(dependency_views.size() == n_dependencies);
    AMANZI_ASSERT(result_views.size() == 1);

    wc = result_views[0];
    pd = dependency_views[0];
    dens = dependency_views[1];
  }

  KeyVector my_keys() const { return { wc_key_ }; }
  KeyVector dependencies() const { return { pd_key_, dens_key_ }; }

  // the model
  KOKKOS_INLINE_FUNCTION void operator()(const int i) const
  {
    wc(i,0) =  pd(i,0) > 0. ? pd(i,0)*dens(i,0) : 0.;
  }

  // derivatives
  //
  // NOTE: the order of these function tags (i.e. Deriv<0>, ...) is set by the
  // above call to dependencies().  Deriv<I> must correspond to the derivative
  // with respect to dependencies()[I];

  // d/dB
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<0>, const int i) const
  {
    //    wc(i,0) = pd(i,0) > 0. ? dens(i,0) : 0.;
    // NOTE this is hacked for surface water stand-alone
    wc(i,0) = dens(i,0);
  }

  KOKKOS_INLINE_FUNCTION void operator()(Deriv<1>, const int i) const
  {
    wc(i,0) = pd(i,0) > 0. ? pd(i,0) : 0.;
  }
  
 private:
  OutView_type wc;
  InView_type pd, dens;

  Key wc_key_;
  Key pd_key_;
  Key dens_key_;

  static Utils::RegisteredFactory<Evaluator, EvaluatorModel_CompositeVector<SurfaceWaterContent>> reg_;

};




} // namespace Relations
} // namespace Flow
} // namespace ATS
