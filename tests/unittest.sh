#!/bin/bash

echo "Unit-testing script for flusi/mhd pseudospectators."
echo

# list all the test scrits you want, separated by spaces
tests=(solid_model.sh jerry.sh insect.sh sphere.sh sphere_restart.sh vortex_ring.sh mhdorszagtang.sh)

# link flusi and mhd from .. to . if this isn't already done.
if [ ! -f flusi ]; then
    ln -s ../flusi .
fi
if [ ! -f mhd ]; then
    ln -s ../mhd .
fi

numtests=0
numsuccess=0
numfail=0

# Get time as a UNIX timestamp (seconds elapsed since Jan 1, 1970 0:00 UTC)
T="$(date +%s)"

for test in ${tests[*]}
do
    numtests=$(($numtests + 1))

    logfile=${test%%.sh}.log
    rm -f $logfile
    touch $logfile

    # run the test, copy the output to the logfile
    ./$test > $logfile
    success=$?
    if [ $success == 0 ]
    then
	echo -e "OK\tRunning test: "${test}", log in: "${logfile}
	numsuccess=$(($numsuccess + 1))
    else
	echo -e "FAIL\tRunning test: "${test}", log in: "${logfile}
	numfail=$(($numfail + 1))
    fi

done

echo
T="$(($(date +%s)-T))"
echo "Time used in tests: ${T} seconds"

echo
echo "Total number of tests: " $numtests
echo "Total number successful: " $numsuccess
echo "Total number failed: " $numfail

if [ $numfail -gt 0 ]
then
    echo "NOT ALL TESTS PASSED: DANGER! DANGER!"
fi
