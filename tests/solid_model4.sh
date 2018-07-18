#!/bin/bash

#-------------------------------------------------------------------------------
# FLUSI (FSI) unit test
# This file contains one specific unit test, and it is called by unittest.sh
#-------------------------------------------------------------------------------
# complete insect test (time consuming but worthwhile)
#-------------------------------------------------------------------------------

# what parameter file
dir=solid_model4_nonconst/
params=$dir"CSM3.ini"

happy=0
sad=0

${mpi_serial} ./flusi --solid ${params}

echo "============================"
echo "run done, analyzing data now"
echo "============================"



#-------------------------------------------------------------------------------
#                               time series
#-------------------------------------------------------------------------------
file="beam_data1.t"

echo comparing $file time series...

${mpi_serial} ./flusi --postprocess --compare-timeseries $file $dir$file

result=$(cat return); rm -f return
if [ "$result" == "0" ]; then
  echo -e ":) Happy, time series: this looks okay! " $file
  happy=$((happy+1))
else
  echo -e ":[ Sad, this is failed! " $file
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
