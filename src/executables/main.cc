/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <iostream>
#include <filesystem>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_YamlParser_decl.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_DefaultComm.hpp"

#include "Epetra_Comm.h"
#include "Epetra_MpiComm.h"
#include "Epetra_SerialComm.h"

#include "VerboseObject_objs.hh"

#include "ats_version.hh"
#include "amanzi_version.hh"
#include "tpl_versions.h"

#include "dbc.hh"
#include "errors.hh"
#include "Key.hh"
#include "ats_driver.hh"

// registration files
#include "ats_registration_files.hh"
#include "PK_Factory.hh"
#include "Evaluator_Factory.hh"

int
main(int argc, char* argv[])
{
#ifdef AMANZI_USE_FENV
  feraiseexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  Teuchos::GlobalMPISession mpiSession(&argc, &argv, 0);
  int rank = mpiSession.getRank();
  Kokkos::initialize();
  int ret = 0;

  {
    std::string input_filename;
    if ((argc >= 2) && (argv[argc - 1][0] != '-')) {
      input_filename = std::string(argv[argc - 1]);
      argc--;
    }

    Teuchos::CommandLineProcessor clp;
    clp.setDocString(
      "Run ATS simulations for ecosystem hydrology.\n\nStandard usage: ats input.xml\n");

    std::string opt_xml_input_filename = "";
    clp.setOption("xml_file", &opt_xml_input_filename, "XML input file");
    std::string opt_input_filename = "";
    clp.setOption("input_file", &opt_input_filename, "input file");

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

    bool list_evals(false);
    clp.setOption(
      "list_evaluators", "no_list_evaluators", &list_evals, "List available evaluators and exit.");

    bool list_pks(false);
    clp.setOption("list_pks", "no_list_pks", &list_pks, "List available PKs and MPCs and exit.");

    clp.throwExceptions(false);
    clp.recogniseAllOptions(true);

    auto parseReturn = clp.parse(argc, argv);
    if (parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
      Kokkos::finalize();
      return 0;
    }
    if (parseReturn != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
      Kokkos::finalize();
      return 1;
    }

#define XSTR(s) STR(s)
#define STR(s) #s
    // check for version info request
    if (version) {
      if (rank == 0) {
        std::cout << "ATS version " << XSTR(ATS_VERSION) << std::endl;
      }
      Kokkos::finalize();
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
      Kokkos::finalize();
      return 0;
    }

    if (list_evals) {
      if (rank == 0) {
        Amanzi::Evaluator_Factory fac;
        std::cout << "Evaluators: "; // no endline is intentional!
        fac.WriteChoices(std::cout);
      }
      Kokkos::finalize();
      return 0;
    }
    if (list_pks) {
      if (rank == 0) {
        Amanzi::PKFactory fac;
        std::cout << "PKs and MPCs: "; // no endline is intentional!
        fac.WriteChoices(std::cout);
      }
      Kokkos::finalize();
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
      Kokkos::finalize();
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
    if (input_filename.empty()) {
      if (!opt_input_filename.empty()) {
        input_filename = opt_input_filename;
      } else if (!opt_xml_input_filename.empty()) {
        input_filename = opt_xml_input_filename;
      } else {
        if (rank == 0) {
          std::cerr << "ERROR: no input file provided" << std::endl;
          clp.printHelpMessage("ats", std::cerr);
        }
        Kokkos::finalize();
        return 1;
      }
    }

    if (!std::filesystem::exists(input_filename)) {
      if (rank == 0) {
        std::cerr << "ERROR: input file \"" << input_filename << "\" does not exist." << std::endl;
      }
      Kokkos::finalize();
      return 1;
    }

    // run the simulation
    // -- create communicator
    auto comm = Amanzi::getDefaultComm();
    auto teuchos_comm = Teuchos::DefaultComm<int>::getComm();

    // -- parse input file
    Teuchos::RCP<Teuchos::ParameterList> plist;
    if (Amanzi::Keys::ends_with(input_filename, ".yaml") ||
        Amanzi::Keys::ends_with(input_filename, ".YAML")) {
      plist = Teuchos::YAMLParameterList::parseYamlFile(input_filename);
    } else {
      plist = Teuchos::getParametersFromXmlFile(input_filename);
    }

    // -- set default verbosity level
    Teuchos::RCP<Teuchos::FancyOStream> fos;
    Teuchos::EVerbosityLevel verbosity_from_list;
    Teuchos::readVerboseObjectSublist(&*plist, &fos, &verbosity_from_list);
    if (verbosity_from_list != Teuchos::VERB_DEFAULT)
      Amanzi::VerboseObject::global_default_level = verbosity_from_list;
    if (!verbosity.empty() ) Amanzi::VerboseObject::global_default_level = opt_level;

    if (Amanzi::VerboseObject::global_default_level != Teuchos::VERB_NONE && (rank == 0)) {
      std::cout
        << "ATS version " << XSTR(ATS_VERSION) << ", Amanzi version " << XSTR(AMANZI_VERSION)
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
        if (rank == 0) {
          std::cerr << "ERROR:" << std::endl << s << std::endl;
        }
        return 1;
      } catch (int& ierr) {
        if (rank == 0) {
          std::cerr << "ERROR: unknown error code " << ierr << std::endl;
        }
        return ierr;
      }
    }
    Teuchos::TimeMonitor::summarize(teuchos_comm.ptr(), std::cout);
  }
  Kokkos::finalize();
  return ret;
}
