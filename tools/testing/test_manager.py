"""
Capability to manage and run ATS system level tests including
regression tests and user demos.

With inspiration, and in some places just stolen, from Ben Andre's
PFloTran regression test suite.

Author: Ethan Coon (ecoon@lanl.gov)
"""
from __future__ import print_function
from __future__ import division
import sys

if sys.hexversion < 0x02070000:
    print(70*"*")
    print("ERROR: The regression test manager requires python >= 2.7.x")
    print("or python 3.")
    print("It appears that you are running python {0}.{1}.{2}".format(
        sys.version_info[0], sys.version_info[1], sys.version_info[2]))
    print(70*"*")
    sys.exit(1)

import datetime
import os
import glob
#import pprint
import shutil
import subprocess
import textwrap
import time
import traceback
import distutils.spawn
import numpy
import collections
import h5py

aliases = dict()


def get_git_hash(directory):
    """Helper function to return the git hash of a repo."""
    import subprocess
    return subprocess.check_output(['git', 'describe', '--always'], cwd=directory).strip()

def version(executable):
    """Helper function to pull a version number from an executable.
    
    For now this is hard-coded and assumes standard ATS locations,
    which is VERY fragile, and may even provide incorrect info.
    This needs to be fixed to take the output of ats --version
    (which doesn't yet work, see ats ticket #12
    """

    try:
        ats_label = get_git_hash(os.environ['ATS_SRC_DIR'])
    except KeyError:
        ats_label = 'unknown'

    try:
        amanzi_label = get_git_hash(os.environ['AMANZI_SRC_DIR'])
    except KeyError:
        amanzi_label = 'unknown'
    v = "ATS: {0}\nAmanzi: {1}".format(ats_label, amanzi_label)
    return v



class NoCatchException(Exception):
    def __init__(self,*args,**kwargs):
        Exception.__init__(self,*args,**kwargs)

if sys.version_info[0] == 2:
    from ConfigParser import SafeConfigParser as config_parser
else:
    from configparser import ConfigParser as config_parser


class TestStatus(object):
    """
    Simple class to hold status info.
    """
    def __init__(self):
        self.fail = 0
        self.warning = 0
        self.error = 0
        self.skipped = 0
        self.test_count = 0

    def __str__(self):
        message = "fail = {0}\n".format(self.fail)
        message += "warning = {0}\n".format(self.warning)
        message += "error = {0}\n".format(self.error)
        message += "skipped = {0}\n".format(self.skipped)
        message += "test_count = {0}\n".format(self.test_count)
        return message


class RegressionTest(object):
    """
    Class to collect data about a test problem, run the problem, and
    compare the results to a known result.
    """
    # things to compare
    _TIME = "time"
    _TIMESTEPS = "timesteps"
    _DEFAULT = "default" # compare a variable
    
    # ways to compare it
    _ABSOLUTE = "absolute"
    _RELATIVE = "relative"
    _PERCENT = "percent"
    _DISCRETE = "discrete"

    # ways to measure what to compare
    _NORM = numpy.inf

    _SUCCESS = 0
    _RESERVED = [_TIME, _TIMESTEPS]
    
    def __init__(self, executable='ats', mpiexec=None, version=None, suffix=None):
        if suffix is not None:
            self._suffix = "."+suffix
        else:
            self._suffix = ".regression"
        
        # misc test parameters
        #        self._pprint = pprint.PrettyPrinter(indent=2)
        self._txtwrap = textwrap.TextWrapper(width=78, subsequent_indent=4*" ")
        self._executable = executable
        self._mpiexec = mpiexec
        self._version = version
        self._input_arg = "--xml_file="
        self._input_suffix = "xml"
        self._np = None
        self._timeout = 162.0
        self._skip_check_gold = False
        self._check_performance = False
        self._num_failed = 0
        self._test_name = None

        # assign default tolerances for different classes of variables
        # absolute min and max thresholds for determining whether to
        # compare to baseline, i.e. if (min_threshold <= abs(value) <=
        # max_threshold) then compare values. By default we use the
        # python definitions for this platform
        self._default_tolerance = {}
        self._default_tolerance[self._TIME] = [1.e-12, self._ABSOLUTE, \
                                               0.0, sys.float_info.max]
        self._default_tolerance[self._DISCRETE] = [0, self._ABSOLUTE, 0, sys.maxsize]
        
        self._default_tolerance[self._DEFAULT] = [1.0e-10, self._RELATIVE, \
                                                  0.0, sys.float_info.max]
        self._eps = 1.e-14

        # dictionary indexed by domain, each of which is a dictionary of
        # (key, tolerance) pairs
        self._criteria = {}
        self._file_prefix = "checkpoint"
        self._file_suffix = ".h5"
        self._file_format = self._file_prefix + "{0:05d}" + self._file_suffix

    def __str__(self):
        message = "  {0} :\n".format(self.name())
        message += "    timeout = {0}\n".format(self._timeout)
        message += "    np = {0}\n".format(self._np)
        message += "    ATS args :\n"
        message += "        exec  : {0}\n".format(self._executable)
        message += "        input : {0}../{1}.{2}\n".format(
            self._input_arg,self._test_name,self._input_suffix)
        for domain,criteria in self._criteria.items():
            message += "    test criteria on \"{0}\" :\n".format(domain)
            for tol in criteria:
                message += "        {}".format(tol)
        return message

    def setup(self, cfg_criteria, test_data, timeout, check_performance, testlog):
        """
        Setup the test object

        cfg_criteria - dict from cfg file, all tests in file
        test_data - dict from cfg file, test specific
        timeout - list(?) from command line option
        check_performance - bool from command line option
        """
        self._test_name = test_data.pop('name')
        self._np = test_data.pop('np', None)
        self._skip_check_gold = test_data.pop('skip_check_gold', None)
        self._check_performance = check_performance

        # timeout : preference (1) command line (2) test data (3) class default
        self._timeout = float(cfg_criteria.pop('timeout', self._timeout))
        self._timeout = float(test_data.pop('timeout', self._timeout))
        if timeout is not None:
            self._timeout = float(timeout[0])

        # pop a default and/or default discrete
        default = test_data.pop(self._DEFAULT,
                                cfg_criteria.pop(self._DEFAULT, None))
        if default is not None:
            tol = self._validate_tolerance(self._DEFAULT, default)
            self._default_tolerance[self._DEFAULT] = tol

        discrete = test_data.pop(self._DISCRETE,
                                cfg_criteria.pop(self._DISCRETE, None))
        if discrete is not None:
            tol = self._validate_tolerance(self._DISCRETE, discrete)
            self._default_tolerance[self._DISCRETE] = tol
            
        # requested criteria, skipping time and timesteps
        for key in set(list(cfg_criteria.keys()) + list(test_data.keys())):
            if key != self._TIME and key != self._TIMESTEPS:
                self._set_criteria(key, cfg_criteria, test_data)

        # always check that write simulation time are the same
        self._set_criteria(self._TIME, cfg_criteria, test_data)

    def name(self):
        return self._test_name

    def dirname(self, gold=False):
        suffix = self._suffix
        if gold:
            suffix = suffix + ".gold"
        return self.name() + suffix        

    def filenames(self, dirname):
        fname_list = sorted(glob.glob(os.path.join(dirname, self._file_prefix+"*"+self._file_suffix)))
        final_name = os.path.join(dirname, 'checkpoint_final.h5')
        if final_name in fname_list:
            fname_list.remove(final_name)
        return fname_list

    def run(self, dry_run, status, testlog):
        """
        Run the test.
        """
        self._cleanup_generated_files()
        self._run_test(self.name(), dry_run, status,
                       testlog)

    def _run_test(self, test_name, dry_run, status, testlog):
        """
        Build up the run command, including mpiexec, np, ats,
        input file, output file. Then run the job as a subprocess.

        * NOTE(bja) - starting in python 3.3, we can use:

          subprocess.Popen(...).wait(timeout)

          to catch hanging jobs, but for python < 3.3 we have to
          manually manage the timeout...?
        """
        command = []
        if self._np is not None:
            if self._mpiexec:
                command.append(self._mpiexec)
                command.append("-np")
                command.append(self._np)
            else:
                # parallel test, but don't have mpiexec, we mark the
                # test as skipped and bail....
                message = self._txtwrap.fill(
                    "WARNING : mpiexec was not provided for a parallel test '{0}'.\n"
                    "This test was skipped!".format(self.name()))
                print(message, file=testlog)
                status.skipped = 1
                return None

        command.append(self._executable)
        command.append("{0}../{1}.{2}".format(self._input_arg,test_name,self._input_suffix))

        test_directory = os.getcwd()
        run_directory = os.path.join(test_directory, self.dirname())
        os.mkdir(run_directory)
        os.chdir(run_directory)
        with open('ats_version.txt', 'w') as fid:
            fid.write(self._version)
        
        print("    cd {0}".format(run_directory), file=testlog)
        print("    {0}".format(" ".join(command)), file=testlog)

        if not dry_run:
            run_stdout = open(test_name + ".stdout", 'w')
            start = time.time()
            proc = subprocess.Popen(command,
                                    shell=False,
                                    stdout=run_stdout,
                                    stderr=run_stdout)
            while proc.poll() is None:
                time.sleep(0.1)
                if time.time() - start > self._timeout:
                    proc.kill()
                    time.sleep(0.1)
                    message = self._txtwrap.fill(
                        "ERROR: job '{0}' has exceeded timeout limit of "
                        "{1} seconds.".format(test_name, self._timeout))
                    print(''.join(['\n', message, '\n']), file=testlog)
            finish = time.time()
            print("    # {0} : run time : {1:.2f} seconds".format(test_name, finish - start), file=testlog)
            ierr_status = abs(proc.returncode)
            run_stdout.close()

            if ierr_status != self._SUCCESS:
                status.fail = 1
                message = self._txtwrap.fill(
                    "FAIL : {name} : {execute} return an error "
                    "code ({status}) indicating the simulation may have "
                    "failed. Please check '{name}.stdout' "
                    "for error messages (included below).".format(
                        execute=self._executable, name=test_name, status=ierr_status))

                print("".join(['\n', message, '\n']), file=testlog)
                print("~~~~~ {0}.stdout ~~~~~".format(test_name), file=testlog)
                try:
                    with open("{0}.stdout".format(test_name), 'r') as tempfile:
                        shutil.copyfileobj(tempfile, testlog)
                except Exception as e:
                    print("   Error opening file: {0}.stdout\n    {1}".format(test_name, e))
                print("~~~~~~~~~~", file=testlog)

        os.chdir(test_directory)

    def _cleanup_generated_files(self):
        """Cleanup old generated files that may be hanging around from a
        previous run.
        """
        run_dir = self.dirname(False)
        shutil.rmtree(run_dir, ignore_errors=True)


    def check(self, status, testlog):
        """
        Check the test results against the gold standard
        """
        self._check_gold(status, testlog)

    def update(self, status, testlog):
        """
        Update the gold standard test results to the current
        output. Both the current regression output and a gold file
        must exist.
        """
        gold_dir = self.dirname(True)
        run_dir = self.dirname(False)
        old_gold_dir = gold_dir+".old"

        # verify that the gold file exists
        if not os.path.isdir(gold_dir):
            print("ERROR: test '{0}' results cannot be updated "
                  "because a gold directory does not "
                  "exist!".format(self.name()), file=testlog)
            status.error = 1
                
        # verify that the regression directory exists
        if not os.path.isdir(run_dir):
            print("ERROR: test '{0}' results cannot be updated "
                  "because no regression run directory "
                  "exists!".format(self.name()), file=testlog)
            status.error = 1

        # check that the regression file was created in the regression directory
        if (not status.error and not status.fail):
            reg_filenames = self.filenames(run_dir)
            if len(reg_filenames) == 0:
                print("ERROR: run for "
                      "test '{0}' did not create the needed regression file(s), "
                      "not updating!".format(self.name(),
                                             reg_name), file=testlog)
                status.fail = 1

        if not status.error and not status.fail:
            # remove old-old gold
            if os.path.isdir(old_gold_dir):
                try:
                    shutil.rmtree(old_gold_dir)
                except Exception as error:
                    message = str(error)
                    message += "\nFAIL : Could not remove old gold '{0}'. ".format(old_gold_dir)
                    message += "Please remove the directory manually!"
                    message += "    rm -rf {0}".format(old_gold_dir)
                    print(message, file=testlog)
                    status.fail = 1

            # move gold -> old_gold
            if not status.fail:
                try:
                    os.rename(gold_dir, old_gold_dir)
                except Exception as error:
                    message = str(error)
                    message += "\nFAIL : Could not back up gold '{0}' as '{1}'. ".format(
                        gold_dir, old_gold_dir)
                    message += "Please rename or remove the gold directory manually!"
                    message += "    mv {0} {1}".format(gold_dir, old_gold_dir)
                    print(message, file=testlog)
                    status.fail = 1

            if not status.fail:
                # move run -> gold
                try:
                    os.rename(run_dir, gold_dir)
                except Exception as error:
                    message = str(error)
                    message += "\nFAIL : Could not move '{0}' to '{1}'. ".format(
                        run_dir, gold_dir)
                    message += "Please rename the directory manually!"
                    message += "    mv {0} {1}".format(run_dir, gold_dir)
                    print(message, file=testlog)
                    status.fail = 1

        print("done", file=testlog)

    def new_test(self, save_dt_history, status, testlog):
        """
        A new test does not have a gold standard regression test. We
        will check to see if a gold standard file exists (an error),
        then create the gold file by copying the current regression
        file to gold.
        """
        gold_dir = self.dirname(True)
        run_dir = self.dirname(False)

        # verify that the gold file does not exist
        if os.path.isdir(gold_dir) or os.path.isfile(gold_dir):
            print("ERROR: test '{0}' was classified as new, "
                               "but a gold file already "
                               "exists!".format(self.name()), file=testlog)
            status.error = 1

        # verify that the run directory exists
        if not os.path.isdir(run_dir):
            print("ERROR: test '{0}' results cannot be new "
                  "because no run directory "
                  "exists!".format(self.name()), file=testlog)
            status.error = 1
        
        # verify that the run directory created good files
        if (not status.error and not status.fail):
            reg_filenames = self.filenames(run_dir)
            if len(reg_filenames) == 0:
                print("ERROR: run for "
                      "test '{0}' did not create the needed regression file "
                      "'{1}', not making new gold!".format(self.name(),
                                                           reg_name), file=testlog)
                status.fail = 1

        if not status.error and not status.fail:
            # move the run to gold
            print("  creating gold directory '{0}'... ".format(self.name()),
                  file=testlog)
            try:
                os.rename(run_dir, gold_dir)
            except Exception as error:
                message = str(error)
                message += "\nFAIL : Could not move '{0}' to '{1}'. ".format(
                    run_dir, gold_dir)
                message += "Please rename the directory manually!"
                message += "    mv {0} {1}".format(run_dir, gold_dir)
                print(message, file=testlog)
                status.fail = 1

            if save_dt_history:
                # this needs some extra things...
                import plot_timestep_history
                import amanzi_xml.utils.io as aio
                import amanzi_xml.utils.search as asearch
                import amanzi_xml.utils.errors as aerrors
                        
                # build the dt history
                with open(os.path.join(gold_dir, self.name()+'.stdout'),'r') as fid:
                    good, bad = plot_timestep_history.parse_logfile(fid)

                if not os.path.isdir('data'):
                    os.mkdir('data')

                with h5py.File(os.path.join('data', "{0}_dts.h5".format(self.name())), 'w') as fid:
                    fid.create_dataset("timesteps", data=86400.*good[:,2]) # note conversion from days back to seconds
                print("Wrote file: data/{0}_dts.h5".format(self.name()), file=testlog)

                # build a new xml file with these specific timestep history
                print("Renaming: {0}.xml to {0}_orig.xml".format(self.name()), file=testlog)
                os.rename(self.name()+".xml", self.name()+"_orig.xml")

                xml = aio.fromFile(self.name()+"_orig.xml", True)

                # -- checkpoint by list of cycles
                filenames = self.filenames(gold_dir)
                chp_cycles = [int(os.path.split(f)[-1][10:-3]) for f in filenames]
                try:
                    chp = asearch.child_by_name(xml, "checkpoint")
                except aerrors.MissingXMLError:
                    raise ValueError("Regression tests are done by checkpointing -- please make sure all tests have at least one checkpoint")

                # empty the checkpoint list
                for i in range(len(chp)):
                    chp.pop(chp[0].get('name'))

                # remove the vis list and coordinator "required times" lists which can affect timestepping
                try:
                    xml.pop("visualization")
                except aerrors.MissingXMLError:
                    pass
                coord = asearch.find_name(xml, "cycle driver")
                try:
                    coord.pop("required times")
                except aerrors.MissingXMLError:
                    pass

                chp.setParameter('cycles', 'Array(int)', chp_cycles)

                # -- update timestep controller, nonlinear solvers
                for ti in asearch.findall_name(xml, "time integrator"):
                    asearch.find_name(ti, "limit iterations").setValue(100)
                    asearch.find_name(ti, "diverged tolerance").setValue(1.e10)

                    asearch.find_name(ti, "timestep controller type").setValue("from file")
                    ts_hist = ti.sublist("timestep controller from file parameters")
                    ts_hist.setParameter("file name", "string", "../data/{0}_dts.h5".format(self.name()))

                # -- write the new xml
                print("Writing: {0}.xml".format(self.name()), file=testlog)
                aio.toFile(xml, '{0}.xml'.format(self.name()))

                # -- run update to get a new run on the dt history
                self.run(False, status, testlog)
                self.update(status, testlog)

        print("done", file=testlog)


    def _check_gold(self, status, testlog):
        """
        Test the output from the run against the known "gold standard"
        output and determine if the test succeeded or failed.

        We return zero on success, one on failure so that the test
        manager can track how many tests succeeded and failed.
        """
        if self._skip_check_gold:
            message = "    Skipping comparison to regression gold file (only test if model runs to completion)."
            print("".join(['\n', message, '\n']), file=testlog)
            return

        # get all gold checkpoint files
        gold_dirname = self.dirname(True)
        gold_files = self.filenames(gold_dirname)

        # get all regression checkpoint files
        reg_dirname = self.dirname(False)
        reg_files = self.filenames(reg_dirname)

        # check that at least one gold file exists
        if len(gold_files) == 0:
            message = self._txtwrap.fill(
                "FAIL: could not find checkpoint files in regression "
                "test gold directory "
                "'{0}'. If this is a new test, please create "
                "it with '--new-test'.".format(gold_dirname))
            print("".join(['\n', message, '\n']), file=testlog)
            status.fail = 1
            return

        # check that gold and regression have same number of files
        if len(gold_files) != len(reg_files):
            message = self._txtwrap.fill(
                "FAIL: differing number of checkpoint files in "
                "gold '{0}' and regression '{1}' directories.".format(gold_dirname, reg_dirname))
            print("".join(['\n', message, '\n']), file=testlog)
            status.fail = 1
            return

        for gold_file, reg_file in zip(gold_files, reg_files):
            try:
                h5_reg = h5py.File(reg_file,'r')
            except Exception as e:
                print("    FAIL: Could not open file: '{0}'".format(reg_file), file=testlog)
                status.fail = 1
                h5_reg = None

            try:
                h5_gold = h5py.File(gold_file,'r')
            except Exception as e:
                print("    FAIL: Could not open file: '{0}'".format(gold_file), file=testlog)
                status.fail = 1
                h5_gold = None

            if h5_reg is not None and h5_gold is not None:
                status.local_fail = 0
                self._compare(h5_reg, h5_gold, status, testlog)
                if status.local_fail == 0:
                    print("    Passed tests {}".format(os.path.split(gold_file)[-1]), file=testlog)
                else:
                    print("    Failed test {}".format(os.path.split(gold_file)[-1]), file=testlog)

            if h5_reg is not None: h5_reg.close()
            if h5_gold is not None: h5_gold.close()
        return
        
    def _compare(self, h5_current, h5_gold, status, testlog):
        """Check that output hdf5 file has not changed from the baseline.
        """
        for key, tolerance in self._criteria.items():
            if key == self._TIME:
                self._check_tolerance(h5_current.attrs['time'], h5_gold.attrs['time'],
                                      self._TIME, tolerance, status, testlog)

            else:
                # find all matches of this key
                reg_matches = []
                if len(key.split('.')) == 1:
                    for k in h5_current.keys():
                        if k.split('.')[0] == key:
                            reg_matches.append(k)
                elif len(key.split('.')) == 2:
                    for k in h5_current.keys():
                        if '.'.join(k.split('.')[0:2]) == key:
                            reg_matches.append(k)
                else:
                    for k in h5_current.keys():
                        if k == key:
                            reg_matches.append(k)

                # find the corresponding matches from gold
                gold_matches = []
                for key in reg_matches:
                    if key in h5_gold.keys():
                        gold_matches.append(key)
                    else:
                        found_match = False
                        key_split = key.split('.')
                        try:
                            my_aliases = aliases[key_split[0]]
                        except KeyError:
                            pass
                        else:
                            for alias in my_aliases:
                                if '.'.join([alias,]+key_split[1:]) in h5_gold.keys():
                                    gold_matches.append(key)
                                    found_match = True
                                    break
                        if not found_match:
                            status.fail = 1
                            print("    FAIL: Cannot find {0} or aliased version in the regression.".format(key),
                                  file=testlog)
                            return

                # we should not get here if this is not true
                assert(len(gold_matches) == len(reg_matches))

                for (g, r) in zip(gold_matches, reg_matches):
                    self._check_tolerance(h5_current[r][:], h5_gold[g][:], r, tolerance, status, testlog)

    def _norm(self, diff):
        """
        Determine the difference between two values
        """
        if type(diff) is numpy.ndarray:
            delta = numpy.linalg.norm(diff.flatten(), self._NORM)
        else:
            delta = abs(diff)
        return delta
                
    def _check_tolerance(self, current, gold, key, tolerance, status, testlog):
        """
        Compare the values using the appropriate tolerance and criteria.
        """
        my_status = 0

        # unpack the tolerance
        tol, tol_type, min_threshold, max_threshold = tuple(tolerance)

        if tol_type == self._ABSOLUTE:
            delta = self._norm(current-gold)

        elif (tol_type == self._RELATIVE or
              tol_type == self._PERCENT):
            if type(gold) is numpy.ndarray:
                min_threshold = max(min_threshold, 1.e-12)
                rel_to = numpy.where(numpy.abs(gold) > min_threshold, gold, min_threshold)
                delta = self._norm((gold - current) / rel_to)

            else:
                delta = self._norm((gold - current) / max(abs(gold), min_threshold))
                
            if tol_type == self._PERCENT:
                delta *= 100.0
        else:
            # should never get here....
            raise RuntimeError("ERROR: unknown test tolerance_type '{0}' for "
                               "variable '{1}, {2}.'".format(tol_type,
                                                          self.name(), key))

        if delta > tol:
            status.fail = 1
            status.local_fail = 1
            print("    FAIL: {0} : {1} : {2} > {3} [{4}]".format(
                    self.name(), key, delta, tol, tol_type), file=testlog)
        else:
            print("    PASS: {0} : {1} : {2} <= {3} [{4}]".format(
                self.name(), key, delta, tol, tol_type), file=testlog)
        return

    def _set_criteria(self, key, cfg_criteria, test_data):
        """
        Our preferred order for selecting test criteria is:
        (1) test data section of the config file
        (2) config-file wide defaults
        (3) hard coded class default
        """
        # parse the criteria
        if key in test_data.keys():
            criteria = key, self._validate_tolerance(key, test_data[key])
        elif key in cfg_criteria.keys():
            criteria = key, self._validate_tolerance(key, cfg_criteria[key])
        elif key == self._TIME:
            criteria = self._TIME, self._default_tolerance[self._TIME]
        else:
            criteria = key, self._default_tolerance[self._DEFAULT]

        if criteria is not None:
            key,tol = criteria
            if tol is not None:
                self._criteria[key] = tol

    def _validate_tolerance(self, key, test_data):
        """
        Validate the tolerance string from a config file.

        Valid input configurations are:
        
        * key = 
        * key = default
        * key = no
        * key = tolerance type
        * key = tolerance type [ min_threshold value [ max_threshold value ] ]

        where the first two use defaults, the third turns the test
        off, and the last two specify the tolerance explicitly, where
        min_threshold and max_threshold are optional

        """
        # deal with defaults first
        if test_data.lower() == "none" or test_data.lower() == "no" or test_data.lower() == "n":
            return None
        if test_data == "" or test_data.lower() == self._DEFAULT:
            return self._default_tolerance[self._DEFAULT]
        if test_data == self._DISCRETE:
            return self._default_tolerance[self._DISCRETE]

        # if we get here, parse the string
        criteria = [None, None, 0.0, sys.float_info.max]
        test_data_s = test_data.split()

        if len(test_data_s) < 2:
            raise RuntimeError("ERROR : Could not convert '{0}' test criteria "
                               "'{1}' into at least a tolerance and type!".format(key, test_data))
        
        try:
            value = float(test_data_s[0])
        except ValueError:
            raise RuntimeError("ERROR : Could not convert '{0}' test criteria "
                               "value '{1}' into a float!".format(key, test_data_s[0]))
        criteria[0] = value

        criteria_type = test_data_s[1]
        if (criteria_type.lower() != self._PERCENT and
            criteria_type.lower() != self._ABSOLUTE and
                criteria_type.lower() != self._RELATIVE):
            raise RuntimeError("ERROR : invalid test criteria string '{0}' "
                               "for '{1}'".format(criteria_type, key))
        criteria[1] = criteria_type

        if (len(test_data_s) > 2):
            try:
                min_threshold = float(test_data_s[2])
            except ValueError:
                raise RuntimeError("ERROR : Could not convert '{0}' test criteria "
                                   "min_threshold value '{1}' into a float!".format(key, test_data_s[2]))
            criteria[2] = min_threshold

        if (len(test_data_s) > 3):
            try:
                max_threshold = float(test_data_s[3])
            except ValueError:
                raise RuntimeError("ERROR : Could not convert '{0}' test criteria "
                                   "max_threshold value '{1}' into a float!".format(key, test_data_s[3]))
            criteria[3] = max_threshold

        return criteria


class RegressionTestManager(object):
    """
    Class to open a configuration file, process it into a group of
    tests, and manage running the tests.

    Notes:

        * The ConfigParser class converts all section names and key
          names into lower case. This means we need to preprocess user
          input names to lower case.

    """

    def __init__(self, executable=None, mpiexec=None, suffix=None):
        self._executable = executable
        self._mpiexec = mpiexec
        self._version = version(executable)
        self._debug = False
        self._file_status = TestStatus()
        self._config_filename = None
        self._default_test_criteria = {}
        self._available_tests = collections.OrderedDict() # ordered to preserve ordering in config files
        self._available_suites = {}
        self._tests = []
        self._txtwrap = textwrap.TextWrapper(width=78, subsequent_indent=4*" ")
        self._suffix = suffix
        #        self._pprint = pprint.PrettyPrinter(indent=2)

    def __str__(self):
        data = "Regression Test Manager :\n"
        data += "    configuration file : {0}\n".format(self._config_filename)
        data += "    default test criteria :\n"
        data += self._dict_to_string(self._default_test_criteria)
        data += "    suites :\n"
        data += self._dict_to_string(self._available_suites)
        data += "    available tests :\n"
        data += self._dict_to_string(self._available_tests)

        data += "Tests :\n"
        for test in self._tests:
            data += test.__str__()

        return data

    def debug(self, debug):
        self._debug = debug

    def num_tests(self):
        return len(self._tests)

    def generate_tests(self, config_file, user_suites, user_tests,
                       timeout, check_performance, testlog):

        """
        Read the config file, validate the input and generate the test objects.
        """
        self._read_config_file(config_file)
        self._validate_suites()
        user_suites, user_tests = self._validate_user_lists(user_suites,
                                                            user_tests, testlog)
        self._create_tests(user_suites, user_tests, timeout, check_performance, testlog)

    def run_tests(self, dry_run, update, new_test, check_only, run_only, testlog, save_dt_history=False):
        """
        Run the tests specified in the config file.

        * dry_run - flag indicates that the test is setup then print
          the command that would be used, but don't actually run
          anything or compare results.

        * new_test - flag indicates that the test is a new test, and
          there should not be a gold standard regression file
          present. Run the executable and create the gold file.

        * update - flag indicates that the output from ats has
          changed, and we want to update the gold standard regression
          file to reflect this. Run the executable and replace the
          gold file.

        * check_only - flag to indicate just diffing the existing
          regression files without rerunning ats.
        """
        if self.num_tests() > 0:
            if new_test:
                self._run_new(save_dt_history, dry_run, testlog)
            elif update:
                self._run_update(dry_run, testlog)
            elif check_only:
                self._check_only(dry_run, testlog)
            elif run_only:
                self._run_only(dry_run, testlog)
            else:
                self._run_check(dry_run, testlog)
        else:
            self._file_status.test_count = 0

    def _run_only(self, dry_run, testlog):
        """
        Run the test and check the results.
        """
        if dry_run:
            print("Dry run:")
        print("Running tests from '{0}':".format(self._config_filename), file=testlog)
        print(50 * '-', file=testlog)

        for test in self._tests:
            status = TestStatus()
            self._test_header(test.name(), testlog)

            test.run(dry_run, status, testlog)

            self._add_to_file_status(status)

            self._test_summary(test.name(), status, dry_run,
                               "passed", "failed", testlog)

        self._print_file_summary(dry_run, "passed", "failed", testlog)

    def _run_check(self, dry_run, testlog):
        """
        Run the test and check the results.
        """
        if dry_run:
            print("Dry run:")
        print("Running tests from '{0}':".format(self._config_filename), file=testlog)
        print(50 * '-', file=testlog)

        for test in self._tests:
            status = TestStatus()
            self._test_header(test.name(), testlog)

            test.run(dry_run, status, testlog)

            if not dry_run and status.skipped == 0:
                test.check(status, testlog)

            self._add_to_file_status(status)

            self._test_summary(test.name(), status, dry_run,
                               "passed", "failed", testlog)

        self._print_file_summary(dry_run, "passed", "failed", testlog)

    def _check_only(self, dry_run, testlog):
        """
        Recheck the regression files from a previous run.
        """
        if dry_run:
            print("Dry run:")
        print("Checking existing test results from '{0}':".format(
            self._config_filename), file=testlog)
        print(50 * '-', file=testlog)

        for test in self._tests:
            status = TestStatus()
            self._test_header(test.name(), testlog)

            if not dry_run and status.skipped == 0:
                test.check(status, testlog)

            self._add_to_file_status(status)

            self._test_summary(test.name(), status, dry_run,
                               "passed", "failed", testlog)

        self._print_file_summary(dry_run, "passed", "failed", testlog)

    def _run_new(self, save_dt_history, dry_run, testlog):
        """
        Run the tests and create new gold files.
        """
        if dry_run:
            print("Dry run:")

        print("New tests from '{0}':".format(self._config_filename), file=testlog)
        print(50 * '-', file=testlog)

        for test in self._tests:
            status = TestStatus()
            self._test_header(test.name(), testlog)
	    
            test.run(dry_run, status, testlog)

            if not dry_run and status.skipped == 0:
                test.new_test(save_dt_history, status, testlog)
            self._add_to_file_status(status)
            self._test_summary(test.name(), status, dry_run,
                               "created", "error creating new test files.", testlog)

        self._print_file_summary(dry_run, "created", "could not be created", testlog)

    def _run_update(self, dry_run, testlog):
        """
        Run the tests and update the gold file with the current output
        """
        if dry_run:
            print("Dry run:")
        print("Updating tests from '{0}':".format(self._config_filename), file=testlog)
        print(50 * '-', file=testlog)

        for test in self._tests:
            status = TestStatus()
            self._test_header(test.name(), testlog)
            test.run(dry_run, status, testlog)

            if not dry_run and status.skipped == 0:
                test.update(status, testlog)
            self._add_to_file_status(status)
            self._test_summary(test.name(), status, dry_run,
                               "updated", "error updating test.", testlog)

        self._print_file_summary(dry_run, "updated", "could not be updated", testlog)

    def _test_header(self, name, testlog):
        """
        Write a header to the log file to separate tests.
        """
        print(40 * '-', file=testlog)
        print("{0}... ".format(name), file=testlog)

    def _test_summary(self, name, status, dry_run,
                      success_message, fail_message, testlog):
        """
        Write the test status information to stdout and the test log.
        """
        if dry_run:
            print("D", end='', file=sys.stdout)
            print(" dry run.", file=testlog)
        else:
            if (status.fail == 0 and
                    status.warning == 0 and
                    status.error == 0 and
                    status.skipped == 0):
                print(".", end='', file=sys.stdout)
                print("{0}... {1}.".format(name, success_message), file=testlog)
            elif status.fail != 0:
                print("F", end='', file=sys.stdout)
                print("{0}... {1}.".format(name, fail_message), file=testlog)
            elif status.warning != 0:
                print("W", end='', file=sys.stdout)
            elif status.error != 0:
                print("E", end='', file=sys.stdout)
            elif status.skipped != 0:
                print("S", end='', file=sys.stdout)
                print("{0}... skipped.".format(name), file=testlog)
            else:
                print("?", end='', file=sys.stdout)

        sys.stdout.flush()

    def _print_file_summary(self, dry_run, success_message, fail_message, testlog):
        """
        Print a summary of the results for this config file
        """
        print("", file=sys.stdout)
        print(50 * '-', file=testlog)
        if dry_run:
            print("{0} : dry run.".format(self._config_filename), file=testlog)
        elif self._file_status.test_count == 0:
            print("{0} : no tests run.".format(self._config_filename), file=testlog)
        else:
            line = "{0} : {1} tests : ".format(self._config_filename,
                                               self._file_status.test_count)
            if self._file_status.fail > 0:
                line = "{0} {1} tests {2}, ".format(
                    line, self._file_status.fail, fail_message)
            if self._file_status.skipped > 0:
                line = "{0} {1} tests {2}, ".format(
                    line, self._file_status.skipped, "skipped")
            num_passed = (self._file_status.test_count -
                          self._file_status.fail - self._file_status.skipped)
            line = "{0} {1} tests {2}".format(line, num_passed, success_message)
            print(line, file=testlog)

    def _add_to_file_status(self, status):
        """
        Add the current test status to the overall status for the file.
        """
        self._file_status.fail += status.fail
        self._file_status.warning += status.warning
        self._file_status.error += status.error
        self._file_status.skipped += status.skipped
        self._file_status.test_count += 1

    def run_status(self):
        return self._file_status

    def display_available_tests(self):
        for test in sorted(self._available_tests.keys()):
            print("{0} {1}".format(os.path.split(os.getcwd())[-1], test))

    def display_available_suites(self):
        for suite in self._available_suites:
            print("    {0} :".format(suite))
            for test in self._available_suites[suite].split():
                print("        {0}".format(test))

    def _read_config_file(self, config_file):
        """
        Read the configuration file.

        Sections : The config file will have known sections:
        "suites", "default-test-criteria".

        All other sections are assumed to be test names.
        """
        if config_file is None:
            raise RuntimeError("Error, must provide a config filename")
        self._config_filename = config_file
        config = config_parser()
        config.read(self._config_filename)

        if config.has_section("default-test-criteria"):
            self._default_test_criteria = \
                self._list_to_dict(config.items("default-test-criteria"))

        if config.has_section("suites"):
            self._available_suites = \
                self._list_to_dict(config.items("suites"))

        self._identify_tests(config)

    def _identify_tests(self, config):
        """
        Create a list of all tests in a config file.

        Assumes every section is a test except for some fixed section
        names

        """
        # section names are test names
        test_names = config.sections()

        # remove the fixed section names
        if config.has_section("default-test-criteria"):
            test_names.remove("default-test-criteria")
        if config.has_section("suites"):
            test_names.remove("suites")

        # all remaining sections should be individual tests
        for test in test_names:
            self._available_tests[test] = self._list_to_dict(config.items(test))
            self._available_tests[test]['name'] = test

    def _dict_to_string(self, data):
        """
        Format dictionary key-value pairs in a string
        """
        temp = ""
        for key, value in data.items():
            temp += "        {0} : {1}\n".format(key, value)
        return temp

    def _list_to_dict(self, input_list):
        """
        Convert a list of key-value pairs into a dictionary.
        """
        output_dict = {}
        for item in input_list:
            output_dict[item[0]] = item[1]
        return output_dict

    def _validate_suites(self):
        """
        Validates the suites defined in configuration file by
        checking that each test in a suite is one of the available
        tests.

        If the config file has an empty suite, we report that to the
        user then remove it from the list.
        """
        invalid_tests = []
        empty_suites = []
        for suite in self._available_suites:
            suite_tests = self._available_suites[suite].split()
            if len(suite_tests) == 0:
                empty_suites.append(suite)
            else:
                # validate the list
                for test in suite_tests:
                    if test not in self._available_tests:
                        name = "suite : '{0}' --> test : '{1}'".format(
                            suite, test)
                        invalid_tests.append(name)

        for suite in empty_suites:
            # empty suite, warn the user and remove it from the list
            del self._available_suites[suite]
            print("DEV WARNING : {0} : cfg validation : empty suite "
                  ": '{1}'".format(self._config_filename, suite))

        if len(invalid_tests) != 0:
            raise RuntimeError("ERROR : suites contain unknown tests in "
                               "configuration file '{0}' : {1}".format(
                                   self._config_filename, invalid_tests))

    def _validate_user_lists(self, user_suites, user_tests, testlog):
        """
        Check that the list of suites or tests passed from the command
        line are valid.
        """
        # if no suites or tests is specified, use all available tests
        if len(user_suites) == 0 and len(user_tests) == 0:
            u_suites = []
            u_tests = self._available_tests
        else:
            # check that the processed user supplied names are valid
            # convert user supplied names to lower case
            u_suites = []
            for suite in user_suites:
                if suite.lower() in self._available_suites:
                    u_suites.append(suite.lower())
                else:
                    message = self._txtwrap.fill(
                        "WARNING : {0} : Skipping requested suite '{1}' (not "
                        "present, misspelled or empty).".format(
                            self._config_filename, suite))
                    print(message, file=testlog)

            u_tests = []
            for test in user_tests:
                if test in self._available_tests:
                    u_tests.append(test.lower())
                else:
                    message = self._txtwrap.fill(
                        "WARNING : {0} : Skipping test '{1}' (not present or "
                        "misspelled).".format(self._config_filename, test))
                    print(message, file=testlog)

        return u_suites, u_tests

    def _create_tests(self, user_suites, user_tests, timeout, check_performance,
                      testlog):
        """
        Create the test objects for all user specified suites and tests.
        """
        all_tests = user_tests
        for suite in user_suites:
            for test in self._available_suites[suite].split():
                all_tests.append(test)

        for test in all_tests:
            #try:
            new_test = RegressionTest(self._executable, self._mpiexec, self._version, self._suffix)
            criteria = self._default_test_criteria.copy()
            new_test.setup(criteria,
                           self._available_tests[test], timeout,
                           check_performance, testlog)
            self._tests.append(new_test)
            #except Exception as error:
            #    raise RuntimeError("ERROR : could not create test '{0}' from "
            #                       "config file '{1}'. {2}".format(
            #                           test, self._config_filename, str(error)))


def config_list_includes_search(options):
    """
    Check if there are any directories in the config list.
    """
    for f in options.configs:
        if (os.path.isdir(f)):
            return True
    return False

def generate_config_file_list(options):
    """
    Try to generate a list of configuration files from the commandline
    options.
    """
    config_file_list = []
    # loop through the list, adding files and searching through directories
    for f in options.configs:
        if not os.path.isabs(f):
            f = os.path.abspath(f)

        if os.path.isfile(f):
            config_file_list.append(f)
        elif os.path.isdir(f):
            search_for_config_files(f, config_file_list)
        else:
            raise RuntimeError("ERROR: specified config file/search "
                               "directory '{0}' does not "
                               "exist!".format(f))

    if options.debug:
        print("\nFound config files:")
        for config_file in config_file_list:
            print("    {0}".format(config_file))

    if len(config_file_list) == 0:
        raise RuntimeError("ERROR: no config files were found. Please specify a "
                           "config file or search directory containing config files.")

    return config_file_list


def search_for_config_files(base_dir, config_file_list):
    """
    recursively search the directory tree, creating a list of all config files
    """
    for root, dirnames, filenames in os.walk(base_dir):
        for filename in filenames:
            if filename.endswith(".cfg") and os.path.isfile(os.path.join(root, filename)):
                config_file_list.append(os.path.join(root, filename))

def check_options(options):
    """
    Run some sanity checks on the commandline options.
    """
    # prevent the user from updating regression output during a
    # recursive search for config files
    if hasattr(options, 'update'):
        if options.update and config_list_includes_search(options):
            raise RuntimeError("ERROR: cannot update gold regression files "
                               "during a recursive search for config files.")

        if options.update and options.new_tests:
            raise RuntimeError("ERROR: cannot create new tests and update gold "
                               "regression files at the same time.")


def check_for_executable(options, testlog):
    """
    Try to verify that we have something reasonable for the executable
    """
    # check the executable
    if options.executable is None:
        # try to detect from env
        try:
            executable = distutils.spawn.find_executable("ats")
        except Exception:
            executable = None
        finally:
            if executable is None:
                options.dry_run = True

    else:
        # absolute path to the executable
        executable = os.path.abspath(options.executable[0])
        # is it a valid file?
        if not os.path.isfile(executable):
            raise RuntimeError("ERROR: executable is not a valid file: "
                               "'{0}'".format(executable))

    if executable is None:
        message = ("\n** WARNING ** : ATS executable was not provided on the command line\n"
                   "               and was not found in the PATH.  Will run as 'dry run.'\n")
        print(message, file=sys.stdout)
        print(message, file=testlog)
        executable = "/usr/bin/false"
    else:
        # try to log some info about executable
        print("ATS information :", file=testlog)
        print("-----------------", file=testlog)
        version_info = version(executable)
        print(version_info, file=testlog)
        print("\n\n", file=testlog)
    return executable


def check_for_mpiexec(options, testlog):
    """
    Try to verify that we have something reasonable for the mpiexec executable
    """

    # check for mpiexec
    mpiexec = None
    if options.mpiexec is None:
        # try to detect from env
        try:
            mpiexec = distutils.spawn.find_executable("mpiexec")
        except IOError:
            mpiexec = None
    else:
        # mpiexec = os.path.abspath(options.mpiexec[0])
        mpiexec = options.mpiexec[0]

    if mpiexec is not None:
        # check that we can use it
        # try to log some info about mpiexec
        print("MPI information :", file=testlog)
        print("-----------------", file=testlog)
        tempfile = "{0}/tmp-ats-regression-test-mpi-info.txt".format(os.getcwd())
        command = [mpiexec, "--version"]
        append_command_to_log(command, testlog, tempfile)
        print("\n\n", file=testlog)
    else:
        message = ("\n** WARNING ** : mpiexec was not provided on the command line.\n"
                   "                All parallel tests will be skipped!\n")
        print(message, file=sys.stdout)
        print(message, file=testlog)

    return mpiexec


def summary_report_by_file(report, outfile):
    """
    Summarize the results for each config file.
    """
    print(70 * '-', file=outfile)
    print("Regression test file summary:", file=outfile)
    for filename in report:
        line = "    {0}... {1} tests : ".format(filename, report[filename].test_count)
        if report[filename].warning > 0:
            line = "{0} {1} test warnings, ".format(line, report[filename].warning)
        if report[filename].error > 0:
            line = "{0} {1} test errors, ".format(line, report[filename].error)

        if report[filename].test_count == 0:
            line = "{0}... no tests were run.".format(line)
        else:
            if report[filename].fail > 0:
                line = "{0} {1} tests failed, ".format(line, report[filename].fail)
            if report[filename].skipped > 0:
                line = "{0} {1} tests skipped, ".format(line, report[filename].skipped)
            if report[filename].fail == 0 and report[filename].skipped == 0:
                line = "{0} all tests passed".format(line)
            else:
                num_passed = (report[filename].test_count - report[filename].fail -
                              report[filename].skipped)
                line = "{0} {1} passed.".format(line, num_passed)

        print("{0}".format(line), file=outfile)

    print("\n", file=outfile)


def summary_report(run_time, report, outfile):
    """
    Overall summary of test results
    """
    print(70 * '-', file=outfile)
    print("Regression test summary:", file=outfile)
    print("    Total run time: {0:4g} [s]".format(run_time), file=outfile)
    test_count = 0
    num_failures = 0
    num_errors = 0
    num_warnings = 0
    num_skipped = 0
    for filename in report:
        test_count += report[filename].test_count
        num_failures += report[filename].fail
        num_errors += report[filename].error
        num_warnings += report[filename].warning
        num_skipped += report[filename].skipped

    print("    Total tests : {0}".format(test_count), file=outfile)

    if num_skipped > 0:
        print("    Skipped : {0}".format(num_skipped), file=outfile)
        success = False

    print("    Tests run : {0}".format(test_count - num_skipped), file=outfile)

    success = True
    if num_failures > 0:
        print("    Failed : {0}".format(num_failures), file=outfile)
        success = False

    if num_errors > 0:
        print("    Errors : {0}".format(num_errors), file=outfile)
        success = False

    if num_warnings > 0:
        print("    Warnings : {0}".format(num_warnings), file=outfile)
        success = False

    if success:
        print("    All tests passed.", file=outfile)

    print("\n", file=outfile)
    return num_failures


def append_command_to_log(command, testlog, tempfile):
    """
    Append the results of a shell command to the test log
    """
    print("$ {0}".format(" ".join(command)), file=testlog)
    testlog.flush()
    with open(tempfile, "w") as tempinfo:
        subprocess.call(command, shell=False,
                        stdout=tempinfo,
                        stderr=subprocess.STDOUT)
        # NOTE(bja) 2013-06 : need a short sleep to ensure the
        # contents get written...?
        time.sleep(0.01)
    with open(tempfile, 'r') as tempinfo:
        shutil.copyfileobj(tempinfo, testlog)
    os.remove(tempfile)    

def setup_testlog(txtwrap, silence=False):
    """
    Create the test log and try to add some useful information about
    the environment.
    """
    now = datetime.datetime.today().strftime("%Y-%m-%d_%H-%M-%S")
    if not os.path.isdir("LOGS"):
        os.mkdir("LOGS")
    filename = os.path.join("LOGS", "ats-tests-{0}.testlog".format(now))

    if not silence:
        print("  Test log file : {0}".format(filename))
    testlog = open(filename, 'w')

    print("ATS Regression Test Log", file=testlog)
    print("Date : {0}".format(now), file=testlog)
    print("System Info :", file=testlog)
    print("    platform : {0}".format(sys.platform), file=testlog)
    test_dir = os.getcwd()
    print("Test directory : ", file=testlog)
    print("    {0}".format(test_dir), file=testlog)

    if 'ATS_SRC_DIR' in os.environ:
        tempfile = "{0}/tmp-ats-regression-test-hg-info.txt".format(test_dir)

        print("\nATS repository status :", file=testlog)
        print("----------------------------", file=testlog)
        if os.path.isdir("{0}/.hg".format(os.environ["ATS_SRC_DIR"])):
            cmd = ["hg", "parent"]
            append_command_to_log(cmd, testlog, tempfile)
            cmd = ["hg", "status", "-q"]
            append_command_to_log(cmd, testlog, tempfile)
            print("\n\n", file=testlog)
        else:
            print("    unknown", file=testlog)

    return testlog

