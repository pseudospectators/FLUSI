#!/bin/bash

#-------------------------------------------------------------------------------
# FLUSI (FSI) unit test
# This file contains one specific unit test, and it is called by unittest.sh
#-------------------------------------------------------------------------------
# SPHERE unit test
# Tests the flow past a sphere at Reynolds number 100
# Using RK2 and vorticity sponge, unlike the second test case ("sphere")
#-------------------------------------------------------------------------------

# set up mpi command (this may be machine dependent!!)
nprocs=$(nproc)
mpi_command="nice -n 19 ionice -c 3 mpiexec --np ${nprocs}"
# what parameter file
params="sphere_sponge_RK2.ini"

happy=0
sad=0

echo "Sphere unit test: phase one"

# list of prefixes the test generates
prefixes=(ux uy uz p vorx vory vorz)
# list of possible times (no need to actually have them)
times=(00000 00100 00200)
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
    reffile=./sphere_sponge_RK2/${p}"_"${t}".ref" 
    
    if [ -f $file ]; then    
        # get four characteristic values describing the field
        ./flusi --postprocess --keyvalues ${file}        
        # and compare them to the ones stored
        if [ -f $reffile ]; then        
            ./flusi --postprocess --compare-keys $keyfile $reffile 
            result=$?
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


#-------------------------------------------------------------------------------
#                               time series
#-------------------------------------------------------------------------------

files=(forces.t)
columns=(1 2 3 4 9 10)
for file in ${files[@]}
do
for col in ${columns[@]}
do
    echo "comparing" $file "column" $col
    new=$(awk -v col=$col '(NR>1) { print  $col }' $file | tail -n1)
    old=$(awk -v col=$col '(NR>1) { print  $col }' sphere_sponge_RK2/$file | tail -n1)
    if [ $new == $old ]; then
      echo ":D HAPPY! timeseries comparison for " $file "column=" $col "succeded"
      happy=$((happy+1))
    else
      echo ":(( Sad: timeseries comparison for " $file " failed"
      sad=$((sad+1))
    fi
done
done


rm -f *.key
rm -f *.h5
rm -f drag_data
rm -f *.t
rm -f runtime*.ini

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
