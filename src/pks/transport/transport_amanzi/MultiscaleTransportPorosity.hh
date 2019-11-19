/*
  Transport PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Base class for continuum multiscale porosity models. 
*/

#ifndef MULTISCALE_TRANSPORT_POROSITY_HH_
#define MULTISCALE_TRANSPORT_POROSITY_HH_

#include <string>
#include <vector>

#include "DenseVector.hh"

namespace Amanzi {
namespace Transport {

class MultiscaleTransportPorosity {
 public:
  virtual ~MultiscaleTransportPorosity() {};

  // Compute solute flux: icomp - component id, phi - matrix porosity,
  // tcc_m_aux - vector of concentration values in secondary nodes,
  // wfm[0|1] - fracture water content at initial and final time moments,
  // wcm[0|1] - matrix water content at initial and final time moments
  virtual double ComputeSoluteFlux(
      double flux_liquid, double& tcc_f, WhetStone::DenseVector& tcc_m, int icomp,
      double dt, double wcf0, double wcf1, double wcm0, double wcm1, double phi) = 0;

  // Number of matrix nodes
  virtual int NumberMatrixNodes() = 0;
};

}  // namespace Transport
}  // namespace Amanzi
  
#endif
  
