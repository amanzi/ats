#!/bin/bash
# usage: collect-failed-logs.sh <tag>
# Run inside the ATS test container. Copy ctest logs and regression-test
# output of the tests that failed in the most recent ctest run into
# /test-logs/<tag>.
shopt -s nullglob
out=/test-logs/$1
mkdir -p "$out"
cd ~/amanzi_builddir/ats
cp Testing/Temporary/LastTest.log "$out"/ 2>/dev/null
[ -f Testing/Temporary/LastTestsFailed.log ] || exit 0
cp Testing/Temporary/LastTestsFailed.log "$out"/

# regression_tests.py is run from the top of ats-regression-tests and writes
# one LOGS/ats-tests-<timestamp>.testlog per invocation (one per ctest test).
# The file name carries no test name, so match on the test name in the content.
logdirs=$(find ~ -type d -path "*/ats-regression-tests/LOGS" 2>/dev/null)

while IFS=: read -r num name; do
  rest=${name#ats_regression_test-}
  [ "$rest" = "$name" ] && continue  # not a regression test; see LastTest.log
  suite=${rest%%-*}
  test=${rest#*-}
  mkdir -p "$out/$suite"
  for logdir in $logdirs; do
    grep -l -F -- "$test" "$logdir"/*.testlog 2>/dev/null | while read -r f; do
      cp "$f" "$out/$suite"/
    done
  done
  # test-local output (stdout captures etc.) lives next to the inputs
  for d in $(find ~ -type d -path "*/ats-regression-tests/$suite" 2>/dev/null); do
    cp "$d/$test"*.stdout "$d/$test"*.regression "$out/$suite"/ 2>/dev/null
  done
done < Testing/Temporary/LastTestsFailed.log
exit 0
