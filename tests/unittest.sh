#!/bin/bash

echo "***************************************************"
echo "Unit-testing script for flusi/mhd pseudospectators."
echo "***************************************************"

# this is to reduce the number of files in the repository
#tar xzf tests_data.tar.gz

# list all the test scrits you want, separated by spaces
tests=( '---insect-module---'
insects_mask_FF_tethered.sh
insects_mask_BB_tethered.sh
insects_mask_FF_tethered_bodyrotated.sh
insects_mask_FF_tethered_kineloader.sh
insects_mask_FF_moving_rotation.sh
insects_mask_FF_moving_rotationtranslation.sh
insects_mask_jerry.sh
insects_mask_corrugation.sh
insects_flow_FF_tethered.sh
insects_flow_BB_tethered_AB2.sh
insects_flow_BB_tethered_AB2_periodic.sh
insects_flow_BB_tethered_RK4.sh
insects_flow_COIN_freeflight.sh
insects_flow_FF_freeflight.sh
'---fluid-without-penalization---'
vortexring1_AB2.sh
vortexring2_RK2.sh
vortexring4_sponge.sh
'---fluid-with-simple-penalization---'
sphere.sh
sphere_restart.sh
sphere_sponge_RK2.sh
'---fluid-structure-interaction---'
solid_model.sh
solid_model2.sh
solid_model3.sh
solid_model4.sh
swimmer1_iteration.sh
swimmer2_staggered.sh
swimmer3_semiimplicit.sh
'---Misc---'
mhdorszagtang.sh
striding.sh
scalar.sh
)

# link flusi and mhd from .. to . if this isn't already done.
if [ ! -f flusi ]; then
    ln -s ../flusi .
fi
if [ ! -f mhd ]; then
    ln -s ../mhd .
fi

if [ -z "$nprocs" ]; then
    nprocs=4
fi
if [ -z "$mpi_command" ]; then
    export mpi_command="nice mpirun -n ${nprocs}"
fi
if [ -z "$mpi_serial" ]; then
    export mpi_serial="nice mpirun -n 1 "
fi

numtests=0
numsuccess=0
numfail=0

echo "employed command for parallel exec: " $mpi_command
echo "employed command for serial   exec: " $mpi_serial
echo "to modify these commands, set \$nprocs, \$mpi_command or \$mpi_serial in shell"

# Get time as a UNIX timestamp (seconds elapsed since Jan 1, 1970 0:00 UTC)
T="$(date +%s)"

for test in ${tests[*]}
do
    numtests=$(($numtests + 1))

    if [[ $test == "---"* ]]; then
    	echo ""
    	echo $test
    	echo ""
  	else
        logfile=${test%%.sh}.log
        rm -f $logfile
        touch $logfile

        # make sure no old data is present before running test
        rm -f *.h5 *.t *.ini

        # run the test, copy the output to the logfile
        ./$test > $logfile
        success=$?
        if [ $success == 0 ]; then
    	     echo -e "OK\tRunning test: "${test}", log in: "${logfile}
    	     numsuccess=$(($numsuccess + 1))
        else
    	     echo -e "FAIL\tRunning test: "${test}", log in: "${logfile}
    	     numfail=$(($numfail + 1))
        fi

        #cleanup
        rm -f *.key *.h5 *.t succ* decom*
        rm -f runtime*
    fi
done

echo
T="$(($(date +%s)-T))"
echo "Time used in tests: ${T} seconds"

echo
echo "Total number of tests: " $(($numsuccess+$numfail))
echo "Total number successful: " $numsuccess
echo "Total number failed: " $numfail

if [ $numfail -gt 0 ]
then
    echo "NOT ALL TESTS PASSED: DANGER! DANGER!"
fi
