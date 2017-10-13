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
subdir=mhdorszagtang
params=${subdir}/mhdorszagtang.ini

happy=0
sad=0

echo "MHD Orszag-Tang case"

# list of prefixes the test generates
prefixes=(ux uy uz bx by bz vorx vory vorz jx jy jz)
# list of possible times (no need to actually have them)
times=(00000 00010)
# run actual test
${mpi_command} ./mhd ${params}

echo "============================"
echo "run done, analyzing data now"
echo "============================"

# loop over all HDF5 files an generate keyvalues using flusi
for p in ${prefixes[@]}
do  
  for t in ${times[@]}
  do
    # *.h5 file coming out of the code
    file=${p}"_"${t}".h5"
    # will be transformed into this *.key file
    keyfile=${p}"_"${t}".key"
    # which we will compare to this *.ref file
    reffile=./${subdir}/${p}"_"${t}".ref" 
    
    if [ -f $file ]; then    
        # get four characteristic values describing the field
        ${mpi_serial} ./flusi --postprocess --keyvalues ${file}        
        # and compare them to the ones stored
        if [ -f $reffile ]; then        
            ${mpi_serial} ./flusi --postprocess --compare-keys $keyfile $reffile 
            result=$(cat return); rm return
            if [ $result == "0" ]; then
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
