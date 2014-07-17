#!/bin/bash

#-------------------------------------------------------------------------------
# FLUSI (FSI) unit test
# This file contains one specific unit test, and it is called by unittest.sh
#-------------------------------------------------------------------------------
# complete insect test (time consuming but worthwhile)
#-------------------------------------------------------------------------------

# what parameter file
params="./solid_model/CSM3.ini"

happy=0
sad=0

./flusi --solid ${params}

echo "============================"
echo "run done, analyzing data now"
echo "============================"



#-------------------------------------------------------------------------------
#                               time series
#-------------------------------------------------------------------------------

files=(beam_data1.t)
columns=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21)
for file in ${files[@]}
do
for col in ${columns[@]}
do
    echo "comparing" $file "column" $col
    new=$(awk -v col=$col '(NR>1) { print  $col }' $file | tail -n1)
    old=$(awk -v col=$col '(NR>1) { print  $col }' solid_model/$file | tail -n1)
    if [ "$new" == "$old" ]; then
      echo ":D HAPPY! timeseries comparison for " $file "column=" $col "succeded"
      happy=$((happy+1))
    else
      echo ":(( Sad: timeseries comparison for " $file " failed"
      sad=$((sad+1))
    fi
    echo " "
done
echo " "
echo " "
echo " "   
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
