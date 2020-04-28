/*
  Copyright 2010-201x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (coonet@ornl.gov)
*/

//! Calculates water depth from pressure.

/*!

.. math:
    h = max\left( \frac{p - p_atm}{\rho g}, 0\right)
  
Parameters include

* `"allow negative`" ``[bool]`` **false** If true, removes the max condition
  and allows negative heights.

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
class SurfaceWaterDepth {
 public:
  static const int n_dependencies = 2;
  static const std::string name;
  
  SurfaceWaterDepth(Teuchos::ParameterList& plist)
  {
    bar_ = plist.get("allow negative in derivative", false);
    
    pd_key_ = Keys::cleanPListName(plist.name());
    std::string domain = Keys::getDomain(pd_key_);
    pres_key_ = Keys::readKey(plist, domain, "pressure", "pressure");
    rho_key_ = Keys::readKey(plist, domain, "mass density", "mass_density_liquid");
  }

  void SetViews(const std::vector<InView_type>& dependency_views,
                const std::vector<OutView_type>& result_views,
                const State& S)
  {
    AMANZI_ASSERT(dependency_views.size() == n_dependencies);
    AMANZI_ASSERT(result_views.size() == 1);

    pd = result_views[0];
    pres = dependency_views[0];
    rho = dependency_views[1];

    p_atm_ = S.Get<double>("atmospheric_pressure", "");
    g_ = S.Get<AmanziGeometry::Point>("gravity", "");
  }

  KeyVector my_keys() const { return { pd_key_ }; }
  KeyVector dependencies() const { return { pres_key_, rho_key_ }; }

  // the model
  KOKKOS_INLINE_FUNCTION void operator()(const int i) const
  {
    double d = (pres(i,0) - p_atm_) / (rho(i,0) * g_[g_.dim()-1]);
    pd(i,0) = d > 0. ? d : 0.;
  }

  // derivatives
  //
  // NOTE: the order of these function tags (i.e. Deriv<0>, ...) is set by the
  // above call to dependencies().  Deriv<I> must correspond to the derivative
  // with respect to dependencies()[I];

  // d/dB
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<0>, const int i) const
  {
    double dddp = 1. / (rho(i,0) * g_[g_.dim()-1]);
    if (bar_) {
      pd(i,0) = dddp;
    } else {
      pd(i,0) = pres(i,0) > p_atm_ ? dddp : 0.;
    }
  }

  // d/dB
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<1>, const int i) const
  {
    double dddp = -(pres(i,0) - p_atm_) * pow((rho(i,0) * g_[g_.dim()-1]), -2);
    if (bar_) {
      pd(i,0) = dddp;
    } else {
      pd(i,0) = pres(i,0) > p_atm_ ? dddp : 0.;
    }      
  }
  
 private:
  OutView_type pd;
  InView_type pres, rho;

  Key pd_key_;
  Key pres_key_;
  Key rho_key_;  

  double p_atm_;
  AmanziGeometry::Point g_;
  bool bar_;
  
  static Utils::RegisteredFactory<Evaluator, EvaluatorModel_CompositeVector<SurfaceWaterDepth>> reg_;

};

} // namespace Relations
} // namespace Flow
} // namespace ATS
