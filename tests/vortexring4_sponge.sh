#!/bin/bash

#-------------------------------------------------------------------------------
# FLUSI (FSI) unit test
# This file contains one specific unit test, and it is called by unittest.sh
#-------------------------------------------------------------------------------
# Vortex ring unit test
# Tests a flow without boundary conditions (no penalization). The initial 
# condition is a vortex ring that travels to the right
#-------------------------------------------------------------------------------

# what parameter file
params="./vortex_ring4_sponge/vortexring4_sponge.ini"

happy=0
sad=0

echo "vortex ring unit test"

# list of prefixes the test generates
prefixes=(ux uy uz p vorx vory vorz)
# list of possible times (no need to actually have them)
times=(000000 000500)
# run actual test
${mpi_command} ./flusi ${params}

echo "============================"
echo "run done, analyzing data now"
echo "============================"

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
    reffile=./vortex_ring4_sponge/${p}"_"${t}".ref" 
    
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
