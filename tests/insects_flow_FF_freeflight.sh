#!/bin/bash

#-------------------------------------------------------------------------------
# FLUSI (FSI) unit test
# This file contains one specific unit test, and it is called by unittest.sh
#-------------------------------------------------------------------------------

# what parameter file

dir="./insects_flow_FF_freeflight/"
params=${dir}"insect_test.ini"
happy=0
sad=0

cp ${dir}/kinematics.ini .

# list of prefixes the test generates
prefixes=(ux uy uz mask vorabs p usx usy usz)
# list of possible times (no need to actually have them)
times=(000000 000250 000500 000750 001000)

# run actual test
${mpi_command} ./flusi ${params}

rm kinematics.ini

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
    reffile=${dir}${p}"_"${t}".ref"

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


#-------------------------------------------------------------------------------
#                               time series
#-------------------------------------------------------------------------------

files=(forces.t forces_part1.t forces_part2.t forces_part3.t kinematics.t energy.t rigidsolidsolver.t ekin.t)

for file in ${files[@]}
do
  echo comparing $file time series...

  ${mpi_serial} ./flusi --postprocess --compare-timeseries $file ${dir}/$file

  result=$(cat return); rm -f return
  if [ "$result" == "0" ]; then
    echo -e ":) Happy, time series: this looks okay! " $file
    happy=$((happy+1))
  else
    echo -e ":[ Sad, time series: this is failed! " $file
    sad=$((sad+1))
  fi
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
