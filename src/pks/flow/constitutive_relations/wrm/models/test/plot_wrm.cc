/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <cstdlib>
#include <iostream>
#include <string>
#include <cmath>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Epetra_SerialComm.h"
#include "Epetra_LocalMap.h"
#include "Epetra_Vector.h"
#include "EpetraExt_VectorOut.h"

#include "wrm_factory.hh"
#include "wrm.hh"


int
main(int argc, char* argv[])
{
  std::cout << "arg count: " << argc << std::endl;
  std::cout << "argv[count-1]: " << argv[argc - 1] << std::endl;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  Teuchos::ParameterList plist;
  std::string xmlFileName = argv[argc - 1];
  std::cout << "reading file: " << xmlFileName << std::endl;
  updateParametersFromXmlFile(xmlFileName, &plist);

  Amanzi::ATS_Physics::Flow::WRMFactory wrmfactory;
  Teuchos::RCP<Amanzi::ATS_Physics::Flow::WRM> wrm = wrmfactory.createWRM(plist);

  // number of fillins between
  int count = 100;

  Epetra_SerialComm comm;
  Epetra_LocalMap map(count, 0, comm);
  Epetra_Vector sat(map);
  Epetra_Vector pc(sat);
  Epetra_Vector krel(sat);

  for (int i = 0; i != count; ++i) {
    sat[i] = i / (count - 1.0);
  }

  double eps = 1e-8;
  bool warned = false;
  for (int i = 0; i != count; ++i) {
    pc[i] = wrm->capillaryPressure(sat[i]);
    if (abs(sat[i] - wrm->saturation(pc[i])) > eps && !warned) {
      std::cout << "ERROR: s != wrm.saturation(wrm.capillaryPressure(s))" << std::endl;
      warned = true;
    }
    krel[i] = wrm->k_relative(pc[i]);
  }

  std::stringstream filename_s_stream;
  filename_s_stream << plist.get<std::string>("wrm type") << "_sat.txt";
  std::string filename_s(filename_s_stream.str());
  for (int j = 0; j != filename_s.length(); ++j)
    if (filename_s[j] == ' ') filename_s[j] = '_';

  std::stringstream filename_pc_stream;
  filename_pc_stream << plist.get<std::string>("wrm type") << "_pc.txt";
  std::string filename_pc(filename_pc_stream.str());
  for (int j = 0; j != filename_pc.length(); ++j)
    if (filename_pc[j] == ' ') filename_pc[j] = '_';

  std::stringstream filename_krel_stream;
  filename_krel_stream << plist.get<std::string>("wrm type") << "_krel.txt";
  std::string filename_krel(filename_krel_stream.str());
  for (int j = 0; j != filename_krel.length(); ++j)
    if (filename_krel[j] == ' ') filename_krel[j] = '_';

  std::cout << "size: " << sat.MyLength() << std::endl;
  EpetraExt::VectorToMatrixMarketFile(filename_s.c_str(), sat, NULL, NULL, false);
  EpetraExt::VectorToMatrixMarketFile(filename_pc.c_str(), pc, NULL, NULL, false);
  EpetraExt::VectorToMatrixMarketFile(filename_krel.c_str(), krel, NULL, NULL, false);
}
