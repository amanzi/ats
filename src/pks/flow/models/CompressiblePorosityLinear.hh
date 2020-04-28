/*
  Copyright 2010-201x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (coonet@ornl.gov)
*/

//! Increase porosity as a function of liquid pressure.

/*!

Evaluates the porosity, given a small compressibility of rock.

Compressible grains are both physically realistic (based on bulk modulus)
and a simple way to provide a non-elliptic, diagonal term for helping
solvers to converge.


Based on a linear increase, i.e.

.. math:
    \phi = \phi_{base} + H(p - p_{atm}) * \alpha

where :math:`H` is the heaviside function and :math:`\alpha` is the provided
compressibility.  If the inflection point is set to zero, the above function
is exact.  However, then the porosity function is not smooth (has
discontinuous derivatives).
  
* `"pore compressibility [Pa^-1]`" ``[double]`` :math:`\alpha` as described above
  
* `"pore compressibility inflection point [Pa]`" ``[double]`` **1000** The inflection point above which the function is linear.

*/

#pragma once

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "Key.hh"
#include "Factory.hh"
#include "StateDefs.hh"
#include "EvaluatorModelByMaterial.hh"
#include "EvaluatorModel_CompositeVector.hh"

namespace ATS {
namespace Flow {
namespace Relations {

using namespace Amanzi;

template <class InView_type, class OutView_type>
class CompressiblePorosityLinear {
 public:
  static const int n_dependencies = 2;
  static const std::string name;

  CompressiblePorosityLinear(Teuchos::ParameterList& plist)
      : alpha_(plist.sublist("model parameters").get<double>("pore compressibility [Pa^-1]")),
        cutoff_(plist.sublist("model parameters").get<double>("pore compressibility inflection point [Pa]", 1000.0))
  {
    phi_key_ = Keys::cleanPListName(plist.name());
    std::string domain = Keys::getDomain(phi_key_);
    phi0_key_ = Keys::readKey(plist, domain, "base porosity", "base_porosity");
    pres_key_ = Keys::readKey(plist, domain, "pressure", "pressure");
  }

  void SetViews(const std::vector<InView_type>& dependency_views,
                const std::vector<OutView_type>& result_views,
                const State& S)
  {
    AMANZI_ASSERT(dependency_views.size() == 2);
    AMANZI_ASSERT(result_views.size() == 1);

    phi = result_views[0];
    phi0 = dependency_views[0];
    pres = dependency_views[1];

    p_atm_ = S.Get<double>("atmospheric pressure", "");
  }

  KeyVector my_keys() const { return { phi_key_ }; }
  KeyVector dependencies() const { return { phi0_key_, pres_key_ }; }

  // the model
  KOKKOS_INLINE_FUNCTION void operator()(const int i) const
  {
    double p_over = pres(i,0) - p_atm_;

    if (p_over > cutoff_) {
      phi(i,0) = phi0(i,0) + alpha_ * ( cutoff_ / 2. + (p_over - cutoff_));
    } else if (p_over > 0.) {
      phi(i,0) = phi0(i,0) + alpha_ * (pow(p_over,2.) / 2. / cutoff_);
    } else {
      phi(i,0) = phi0(i,0);
    }

    phi(i,0) = phi(i,0) > 1.0 ? 1.0 : phi(i,0);
  }

  // derivatives
  //
  // NOTE: the order of these function tags (i.e. Deriv<0>, ...) is set by the
  // above call to dependencies().  Deriv<I> must correspond to the derivative
  // with respect to dependencies()[I];

  // d/d phi0
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<0>, const int i) const
  {
    operator()(i);
    phi(i,0) = phi(i,0) == 1.0 ? 0. : 1.0;
  }

  // d/d pres
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<1>, const int i) const
  {
    operator()(i);
    if (phi(i,0) == 1.0) {
      phi(i,0) = 0.;
    } else {
      double p_over = pres(i,0) - p_atm_;

      if (p_over > cutoff_) {
        phi(i,0) = alpha_;
      } else if (p_over > 0.) {
        phi(i,0) = alpha_ * p_over / cutoff_;
      } else {
        phi(i,0) = 0.;
      }
    }
  }
  
 private:
  OutView_type phi;
  InView_type phi0, pres;

  Key phi_key_;
  Key phi0_key_;
  Key pres_key_;

  double alpha_, cutoff_, p_atm_;

  static Utils::RegisteredFactory<Evaluator, EvaluatorModelByMaterial<CompressiblePorosityLinear>> by_material_reg_;
  static Utils::RegisteredFactory<Evaluator, EvaluatorModel_CompositeVector<CompressiblePorosityLinear>> global_reg_;

};


} // namespace Relations
} // namespace Flow
} // namespace ATS
