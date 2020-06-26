/*
  Copyright 2010-201x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (coonet@ornl.gov)
*/

//! Manning's conductivity provides the nonlinear coefficient for overland flow.

/*!

Note we include the density term here, simply to avoid having to multiply it in
another evaluator.
  
.. math:
    k = n_l h \frac{h^{\beta}{n_{mann} max(sqrt(|S|), \epsilon)}
  
Parameters include:

* `"Manning exponent [-]`" ``[double]`` **2/3** :math:`\beta` above, the exponent.
* `"minimum slope magnitude [-]`" ``[double]`` **0.01** :math:`\epsilon` above,
  keeps the denominator from being zero.

*/

#pragma once

#include <cmath>

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
class ManningConductivity {
 public:
  static const int n_dependencies = 4;
  static const std::string name;
  
  ManningConductivity(Teuchos::ParameterList& plist)
  {
    // note +1 comes from additional leading term h
    beta_ = plist.sublist("model parameters").get("Manning exponent [-]", 2./3) + 1.0;
    eps_ = plist.sublist("model parameters").get("minimum slope magnitude [-]", 0.01);
    
    cond_key_ = Keys::cleanPListName(plist.name());
    std::string domain = Keys::getDomain(cond_key_);
    h_key_ = Keys::readKey(plist, domain, "water depth", "ponded_depth");
    dens_key_ = Keys::readKey(plist, domain, "density", "molar_density_liquid");
    mann_key_ = Keys::readKey(plist, domain, "Manning coeficient", "manning_coefficient");
    slope_key_ = Keys::readKey(plist, domain, "slope magnitude", "slope_magnitude");
  }

  void SetViews(const std::vector<InView_type>& dependency_views,
                const std::vector<OutView_type>& result_views,
                const State& S)
  {
    AMANZI_ASSERT(dependency_views.size() == n_dependencies);
    AMANZI_ASSERT(result_views.size() == 1);

    cond = result_views[0];
    h = dependency_views[0];
    dens = dependency_views[1];
    mann = dependency_views[2];
    slope = dependency_views[3];
  }

  KeyVector my_keys() const { return { cond_key_ }; }
  KeyVector dependencies() const { return { h_key_, dens_key_, mann_key_, slope_key_ }; }

  // the model
  KOKKOS_INLINE_FUNCTION void operator()(const int i) const
  {
    double h0 = h(i) > 0. ? h(i) : 0.;
    double sm = sqrt( slope(i) < eps_ ? eps_ : slope(i) );
    cond(i) = dens(i) * pow(h0, beta_) / mann(i) / sm;
  }

  // derivatives
  //
  // NOTE: the order of these function tags (i.e. Deriv<0>, ...) is set by the
  // above call to dependencies().  Deriv<I> must correspond to the derivative
  // with respect to dependencies()[I];

  KOKKOS_INLINE_FUNCTION void operator()(Deriv<0>, const int i) const
  {
    double h0 = h(i) > 0. ? h(i) : 0.;
    double sm = sqrt( slope(i) < eps_ ? eps_ : slope(i) );
    cond(i) = beta_ * dens(i) * pow(h0, beta_-1.0) / mann(i) / sm;
  }
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<1>, const int i) const
  {
    double h0 = h(i) > 0. ? h(i) : 0.;
    double sm = sqrt( slope(i) < eps_ ? eps_ : slope(i) );
    cond(i) = pow(h0, beta_) / mann(i) / sm;
  }
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<2>, const int i) const
  {
    double h0 = h(i) > 0. ? h(i) : 0.;
    double sm = sqrt( slope(i) < eps_ ? eps_ : slope(i) );
    cond(i) = - dens(i) * pow(h0, beta_) * pow(mann(i),-2) / sm;
  }
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<3>, const int i) const
  {
    double h0 = h(i) > 0. ? h(i) : 0.;
    if (slope(i) < eps_) {
      cond(i) = 0.;
    } else {
      cond(i) = -0.5 * dens(i) * pow(h0, beta_) * pow(slope(i), -1.5) / mann(i);
    }
  }

 private:
  OutView_type cond;
  InView_type dens, h, slope, mann;

  Key cond_key_;
  Key dens_key_;
  Key h_key_;  
  Key slope_key_;  
  Key mann_key_;  

  double beta_;
  double eps_;
  
  static Utils::RegisteredFactory<Evaluator, EvaluatorModel_CompositeVector<ManningConductivity>> reg_;

};

} // namespace Relations
} // namespace Flow
} // namespace ATS
