/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Reads porosity from checkpoint file.

  Authors: Ahmad Jan (jana@ornl.gov)
*/

#include "porosity_fromcheckpointfile_evaluator.hh"

namespace Amanzi {
namespace Flow {

// registry of method 
  Utils::RegisteredFactory<FieldEvaluator,PorosityFromCheckpointFileEvaluator> PorosityFromCheckpointFileEvaluator::fac_("porosity from checkpoint file"); 

} //namespace
} //namespace
