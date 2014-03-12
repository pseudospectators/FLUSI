#!/bin/bash

#-------------------------------------------------------------------------------
# FLUSI (FSI) unit test
# This file contains one specific unit test, and it is called by unittest.sh
#-------------------------------------------------------------------------------
# complete insect test (time consuming but worthwhile)
#-------------------------------------------------------------------------------

# set up mpi command (this may be machine dependent!!)
mpi_command="mpiexec --np 4 -host localhost"
# what parameter file
params="insect.ini"

happy=0
sad=0

echo "big insect test"

# list of prefixes the test generates
prefixes=(mask p usx usy usz ux uy uz vorx vory vorz)
# list of possible times (no need to actually have them)
times=(00000 00010 00020 00030 00040 00050)
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
    # *.h5 file coming out of the code
    file=${p}"_"${t}".h5"
    # will be transformed into this *.key file
    keyfile=${p}"_"${t}".key"
    # which we will compare to this *.ref file
    reffile=./insect/${p}"_"${t}".ref" 
    
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

files=(forces.t forces_part1.t forces_part2.t)
columns=(1 2 4 5 6 7 9 12)
for file in ${files[@]}
do
for col in ${columns[@]}
do
    echo "comparing" $file "column" $col
    new=$(awk -v col=$col '(NR>1) { print  $col }' $file | tail -n1)
    old=$(awk -v col=$col '(NR>1) { print  $col }' insect/$file | tail -n1)
    if [ $new == $old ]; then
      echo ":D HAPPY! timeseries comparison for " $file "column=" $col "succeded"
      happy=$((happy+1))
    else
      echo ":(( Sad: timeseries comparison for " $file " failed"
      sad=$((sad+1))
    fi
done
done


columns=(1 2 4 5 6 7 8 9 12 13 14)
file=kinematics.t
for col in ${columns[@]}
do
    echo "comparing" $file "column" $col
    new=$(awk -v col=$col '(NR>1) { print  $col }' $file | tail -n1)
    old=$(awk -v col=$col '(NR>1) { print  $col }' insect/$file | tail -n1)
    if [ $new == $old ]; then
      echo ":D HAPPY! timeseries comparison for " $file "column=" $col "succeded"
      happy=$((happy+1))
    else
      echo ":(( Sad: timeseries comparison for " $file " failed"
      sad=$((sad+1))
    fi
done


#-------------------------------------------------------------------------------
#                               cleanup
#-------------------------------------------------------------------------------
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