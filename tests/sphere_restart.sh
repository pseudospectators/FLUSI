#!/bin/bash

#-------------------------------------------------------------------------------
# FLUSI (FSI) unit test
# This file contains one specific unit test, and it is called by unittest.sh
#-------------------------------------------------------------------------------
# SPHERE unit test
# Tests the flow past a sphere at Reynolds number 100
# Test has two stages:
#       - flow past a sphere
#       - same flow, but interrupted and restarted
#-------------------------------------------------------------------------------


happy=0
sad=0

#-------------------------------------------------------------------------------
#                               STAGE TWO: restarting test
#-------------------------------------------------------------------------------
echo "Sphere unit test: restart facility test"

# list of prefixes the test generates
prefixes=(ux uy uz p vorx vory vorz)
# list of possible times (no need to actually have them)
times=(000000 002000)
# run first part: starting (runs untill T=1.0)
${mpi_command} ./flusi ./sphere/testing_sphere_start.ini
echo -e "\t\t============================"
echo -e "\t\t RESTARTING"
echo -e "\t\t============================"
# run second part: retake the simulation and run to T=2.0
${mpi_command} ./flusi ./sphere/testing_sphere_restart.ini

# loop over all HDF5 files an generate keyvalues using flusi
for p in ${prefixes[@]}
do
  for t in ${times[@]}
  do
    echo "--------------------------------------------------------------------"
    # *.h5 file coming out of the code
    file=${p}"_"${t}".h5"
    # will be transformed into this *.key file
    keyfile=${p}"_"${t}".key"
    # which we will compare to this *.ref file
    reffile=./sphere/${p}"_"${t}".ref"

    if [ -f $file ]; then
        # get four characteristic values describing the field
        ${mpi_serial} ./flusi --postprocess --keyvalues ${file}
        # and compare them to the ones stored
        if [ -f $reffile ]; then
            ${mpi_serial} ./flusi --postprocess --compare-keys $keyfile $reffile
            result=$(cat return); rm -f return
            if [ "$result" == "0" ]; then
              echo -e ":) Happy, this looks okay! " $keyfile $reffile
              happy=$((happy+1))
            else
              echo -e ":[ Sad, this is failed! " $keyfile $reffile
              sad=$((sad+1))
            fi
        else
            sad=$((sad+1))
            echo -e ":[ Sad: Reference file not found"
        fi
    else
        sad=$((sad+1))
        echo -e ":[ Sad: output file not found"
    fi
    echo " "
    echo " "

  done
done


echo -e "\thappy tests: \t" $happy
echo -e "\tsad tests: \t" $sad

#-------------------------------------------------------------------------------
#                               RETURN
#-------------------------------------------------------------------------------
if [ $sad == 0 ]
then
  exit 0
else
  exit 999
fi
