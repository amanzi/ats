/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluates the porosity, given a small compressibility of rock.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "compressible_porosity_leijnse_evaluator.hh"

namespace Amanzi {
namespace Flow {

// registry of method 
  Utils::RegisteredFactory<Evaluator,CompressiblePorosityLeijnseEvaluator> CompressiblePorosityLeijnseEvaluator::fac_("compressible porosity leijnse"); 

} //namespace
} //namespace
