/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <iostream>

#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include "Epetra_SerialComm.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_DefaultComm.hpp"

#include "VerboseObject_objs.hh"

#include "ats_version.hh"
#include "amanzi_version.hh"
#include "tpl_versions.h"

#include "dbc.hh"
#include "errors.hh"
#include "ats_driver.hh"

// registration files
#include "ats_registration_files.hh"


// include fenv if it exists
#include "boost/version.hpp"
#if (BOOST_VERSION / 100 % 1000 >= 46)
#  include "boost/config.hpp"
#  ifndef BOOST_NO_FENV_H
#    ifdef _GNU_SOURCE
#      define AMANZI_USE_FENV
#      include "boost/detail/fenv.hpp"
#    endif
#  endif
#endif

#include "boost/filesystem.hpp"

int
main(int argc, char* argv[])
{
#ifdef AMANZI_USE_FENV
  //  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  feraiseexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  Teuchos::GlobalMPISession mpiSession(&argc, &argv, 0);
  int rank = mpiSession.getRank();

  std::string input_filename;
  if ((argc >= 2) && (argv[argc - 1][0] != '-')) {
    input_filename = std::string(argv[argc - 1]);
    argc--;
  }

  Teuchos::CommandLineProcessor clp;
  clp.setDocString(
    "Run ATS simulations for ecosystem hydrology.\n\nStandard usage: ats input.xml\n");

  std::string opt_input_filename = "";
  clp.setOption("xml_file", &opt_input_filename, "XML input file");

  bool version(false);
  clp.setOption("version", "no_version", &version, "Print version number and exit.");

  bool print_version(false);
  clp.setOption(
    "print_version", "no_print_version", &print_version, "Print full version info and exit.");

  std::string verbosity;
  clp.setOption("verbosity",
                &verbosity,
                "Default verbosity level: \"none\", \"low\", \"medium\", \"high\", \"extreme\".");

  std::string writing_rank;
  clp.setOption("write_on_rank", &writing_rank, "Rank on which to write VerboseObjects");

  clp.throwExceptions(false);
  clp.recogniseAllOptions(true);

  auto parseReturn = clp.parse(argc, argv);
  if (parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) { return 0; }
  if (parseReturn != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) { return 1; }

#define XSTR(s) STR(s)
#define STR(s) #s
  // check for version info request
  if (version) {
    if (rank == 0) { std::cout << "ATS version " << XSTR(ATS_VERSION) << std::endl; }
    return 0;
  }
  if (print_version) {
    if (rank == 0) {
      std::cout << "ATS version     " << XSTR(ATS_VERSION) << std::endl;
      std::cout << "GIT branch      " << XSTR(ATS_GIT_BRANCH) << std::endl;
      std::cout << "GIT global hash " << XSTR(ATS_GIT_GLOBAL_HASH) << std::endl;
      std::cout << std::endl;
      std::cout << "Amanzi version  " << XSTR(AMANZI_VERSION) << std::endl;
      std::cout << "GIT branch      " << XSTR(AMANZI_GIT_BRANCH) << std::endl;
      std::cout << "GIT global hash " << XSTR(AMANZI_GIT_GLOBAL_HASH) << std::endl;
      std::cout << std::endl;
    }
    return 0;
  }

  // parse the verbosity level
  Teuchos::EVerbosityLevel opt_level;
  if (verbosity.empty()) {
    // pass
  } else if (verbosity == "none") {
    opt_level = Teuchos::VERB_NONE;
  } else if (verbosity == "low") {
    opt_level = Teuchos::VERB_LOW;
  } else if (verbosity == "medium") {
    opt_level = Teuchos::VERB_MEDIUM;
  } else if (verbosity == "high") {
    opt_level = Teuchos::VERB_HIGH;
  } else if (verbosity == "extreme") {
    opt_level = Teuchos::VERB_EXTREME;
  } else {
    if (rank == 0) {
      std::cerr << "ERROR: invalid verbosity level \"" << verbosity << "\"" << std::endl;
      clp.printHelpMessage("ats", std::cerr);
    }
    return 1;
  }

  // parse the writing rank
  if (writing_rank.empty()) {
    // pass
  } else {
    int writing_rank_j;
    try {
      writing_rank_j = std::stoi(writing_rank);
    } catch (std::invalid_argument& e) {
      std::cerr << "ERROR: invalid writing rank \"" << writing_rank << "\"" << std::endl;
      clp.printHelpMessage("ats", std::cerr);
    }
    if (writing_rank_j < 0) {
      std::cerr << "ERROR: invalid writing rank \"" << writing_rank << "\"" << std::endl;
      clp.printHelpMessage("ats", std::cerr);
    }
    Amanzi::VerboseObject::global_writing_rank = writing_rank_j;
  }

  // parse the input file and check validity
  if (input_filename.empty() && !opt_input_filename.empty()) input_filename = opt_input_filename;
  if (input_filename.empty()) {
    if (rank == 0) {
      std::cerr << "ERROR: no input file provided" << std::endl;
      clp.printHelpMessage("ats", std::cerr);
    }
    return 1;
  } else if (!boost::filesystem::exists(input_filename)) {
    if (rank == 0) {
      std::cerr << "ERROR: input file \"" << input_filename << "\" does not exist." << std::endl;
    }
    return 1;
  }

  // run the simulation
  // -- create communicator
  auto comm = Amanzi::getDefaultComm();
  auto teuchos_comm = Teuchos::DefaultComm<int>::getComm();

  // -- parse input file
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(input_filename);

  // -- set default verbosity level
  Teuchos::RCP<Teuchos::FancyOStream> fos;
  Teuchos::EVerbosityLevel verbosity_from_list;
  Teuchos::readVerboseObjectSublist(&*plist, &fos, &verbosity_from_list);
  if (verbosity_from_list != Teuchos::VERB_DEFAULT)
    Amanzi::VerboseObject::global_default_level = verbosity_from_list;
  if (!verbosity.empty()) Amanzi::VerboseObject::global_default_level = opt_level;

  if (Amanzi::VerboseObject::global_default_level != Teuchos::VERB_NONE && (rank == 0)) {
    std::cout << "ATS version " << XSTR(ATS_VERSION) << ", Amanzi version " << XSTR(AMANZI_VERSION)
              << std::endl
              << "================================================================================="
                 "====================="
              << std::endl
              << std::flush;
  }


  // create the top level driver and run simulation
  int ret = 0;
  {
    auto wallclock_timer = Teuchos::TimeMonitor::getNewCounter("wallclock duration");
    ATS::ATSDriver driver(plist, wallclock_timer, teuchos_comm, comm);
    try {
      ret = driver.run();
    } catch (std::string& s) {
      if (rank == 0) { std::cerr << "ERROR:" << std::endl << s << std::endl; }
      return 1;
    } catch (int& ierr) {
      if (rank == 0) { std::cerr << "ERROR: unknown error code " << ierr << std::endl; }
      return ierr;
    }
  }
  Teuchos::TimeMonitor::summarize(teuchos_comm.ptr(), std::cout);
  return ret;
}
