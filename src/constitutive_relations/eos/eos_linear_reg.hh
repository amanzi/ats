/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  ATS

  Linear density/viscosity EOS, defaults to reasonable values for water.

  http://software.lanl.gov/ats/trac

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "eos_linear.hh"

namespace ATS {
namespace Relations {

template <>
const Utils::RegisteredFactory<Evaluator,EvaluatorModel_CompositeVector<EOSLinear<MultiVectorView_type, cMultiVectorView_type> > fac_cv_("eos linear");


// registry of method
Utils::RegisteredFactory<EOS,EOSLinear> EOSLinear::factory_("linear");

} // namespace
} // namespace
