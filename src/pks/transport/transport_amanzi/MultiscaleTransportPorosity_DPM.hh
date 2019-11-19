/*
  Transport PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  A two-scale porosity model (fracture + matrix) aka dual porosity
  model. Current naming convention is that the fields used in the 
  single-porosity model correspond now to the fracture continuum.
  Example: tcc = total component concentration in the fracture continuum;
           tcc_matrix = total component concentration in the matrix continuum.
*/

#ifndef MULTISCALE_TRANSPORT_POROSITY_DPM_HH_
#define MULTISCALE_TRANSPORT_POROSITY_DPM_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "DenseVector.hh"

#include "MultiscaleTransportPorosity.hh"

namespace Amanzi {
namespace Transport {

class MultiscaleTransportPorosity_DPM : public MultiscaleTransportPorosity {
 public:
  MultiscaleTransportPorosity_DPM(const Teuchos::ParameterList& plist);
  ~MultiscaleTransportPorosity_DPM() {};

  // Compute solute flux: icomp - component id, phi - matrix porosity,
  // tcc_m_aux - vector of concentration values in secondary nodes,
  // wfm[0|1] - fracture water content at initial and final time moments,
  // wcm[0|1] - water content at initial and final time moments
  virtual double ComputeSoluteFlux(
      double flux_liquid, double& tcc_f, WhetStone::DenseVector& tcc_m, int icomp,
      double dt, double wcf0, double wcf1, double wcm0, double wcm1, double phi) override;

  // Number of matrix nodes
  virtual int NumberMatrixNodes() override { return 1; }

 private:
  std::vector<double> mol_diff_;
  double warren_root_, tau_, depth_;

  static Utils::RegisteredFactory<MultiscaleTransportPorosity, MultiscaleTransportPorosity_DPM> factory_;
};

}  // namespace Transport
}  // namespace Amanzi
  
#endif
  
