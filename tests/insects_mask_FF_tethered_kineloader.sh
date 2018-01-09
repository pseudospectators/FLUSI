#!/bin/bash

#-------------------------------------------------------------------------------
# FLUSI (FSI) unit test
# This file contains one specific unit test, and it is called by unittest.sh
#-------------------------------------------------------------------------------
# jerry mask test
#-------------------------------------------------------------------------------

# what parameter file
dir="insects_mask_FF_tethered_kineloader/"
params=${dir}"insect_test.ini"
cp ${dir}data_kin.dat ./
happy=0
sad=0


# list of prefixes the test generates
prefixes=(mask usx usy usz)
# list of possible times (no need to actually have them)
times=(000000 000750 001500 002250 003000)
# run actual test
${mpi_command} ./flusi --dry-run ${params}
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
    reffile=./${dir}${p}"_"${t}".ref"

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

rm data_kin.dat


#-------------------------------------------------------------------------------
#                               time series
#-------------------------------------------------------------------------------

file=kinematics.dry-run.t

  echo comparing $file time series...

  ${mpi_serial} ./flusi --postprocess --compare-timeseries $file $dir/$file

  result=$(cat return); rm -f return
  if [ "$result" == "0" ]; then
    echo -e ":) Happy, time series: this looks okay! " $file
    happy=$((happy+1))
  else
    echo -e ":[ Sad, time series: this is failed! " $file
    sad=$((sad+1))
  fi




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
