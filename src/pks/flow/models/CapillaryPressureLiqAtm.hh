/*
  Copyright 2010-201x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (coonet@ornl.gov)
*/

//! Simple capillary pressure, p_atm - p

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
class CapillaryPressureLiqAtm {
 public:
  static const int n_dependencies = 1;
  static const std::string name;

  CapillaryPressureLiqAtm(Teuchos::ParameterList& plist)
  {
    pc_key_ = Keys::cleanPListName(plist.name());
    std::string domain = Keys::getDomain(pc_key_);
    pres_key_ = Keys::readKey(plist, domain, "pressure", "pressure");
  }

  void SetViews(const std::vector<InView_type>& dependency_views,
                const std::vector<OutView_type>& result_views,
                const State& S)
  {
    AMANZI_ASSERT(dependency_views.size() == n_dependencies);
    AMANZI_ASSERT(result_views.size() == 1);

    pc = result_views[0];
    pres = dependency_views[0];

    p_atm_ = S.Get<double>("atmospheric pressure", "");
  }

  KeyVector my_keys() const { return { pc_key_ }; }
  KeyVector dependencies() const { return { pres_key_ }; }

  // the model
  KOKKOS_INLINE_FUNCTION void operator()(const int i) const
  {
    pc(i) =  p_atm_ - pres(i);
  }

  // derivatives
  //
  // NOTE: the order of these function tags (i.e. Deriv<0>, ...) is set by the
  // above call to dependencies().  Deriv<I> must correspond to the derivative
  // with respect to dependencies()[I];

  // d/dB
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<0>, const int i) const
  {
    pc(i) = -1.0;
  }

 private:
  OutView_type pc;
  InView_type pres;

  Key pc_key_;
  Key pres_key_;

  double p_atm_;

  static Utils::RegisteredFactory<Evaluator, EvaluatorModel_CompositeVector<CapillaryPressureLiqAtm>> reg_;

};




} // namespace Relations
} // namespace Flow
} // namespace ATS
