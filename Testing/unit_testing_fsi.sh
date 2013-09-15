#!/bin/bash

if [ ! -f flusi ]; then
  # if no flusi excecutable present, try copying it from root folder
  if [ -f ../flusi ]; then
      cp ../flusi ./
  else
      echo "no flusi excecutable present..."
      exit 1
  fi
fi

# meaningful name for the test
names[1]="Vortex Ring 1: serial, no penalization"
names[2]="Vortex Ring 1: parallel, no penalization"
names[3]="Sphere, with penalization, longtime"
names[4]="restart facility test"

# which params file?
params_files[1]="Testing_VortexRing.ini"
params_files[2]="Testing_VortexRing.ini"
params_files[3]="Testing_Sphere.ini"
params_files[4]="Testing_Sphere.ini"

# mpi or not?
parallel[1]="serial"
parallel[2]="parallel"
parallel[3]="parallel"
parallel[4]="restart"

# where is the reference data?
dirs[1]="VortexRing"
dirs[2]="VortexRing"
dirs[3]="Sphere"
dirs[4]="Sphere"

# how many tests?
n_tests=${#names[@]} 

# list of possible prefixes (no worries if some do not exist for 
# that particular test)
prefixes=(ux uy uz p vorx vory vorz mask usx usy usz)
# list of possible times (no need to actually have them)
times=(00000 00010 00100 00200)

passed=0
error=0
warn=0
checked=0

# Reset
Color_Off='\e[0m'       # Text Reset
Black='\e[0;30m'        # Black
Red='\e[0;31m'          # Red
Green='\e[0;32m'        # Green
Yellow='\e[0;33m'       # Yellow
Blue='\e[0;34m'         # Blue
Purple='\e[0;35m'       # Purple
Cyan='\e[0;36m'         # Cyan
White='\e[0;37m'        # White


for (( i=1; i<=n_tests; i++ ))
do
  echo -e ${Purple} "-----------------------------" ${Color_Off}
  echo -e ${Purple} "Testing: " ${names[i]} " test "${i}"/"${n_tests} ${Color_Off}
  echo -e ${Purple} "-----------------------------" ${Color_Off}
  sleep 0.5
  
  #-----------------------------------------------
  # run actual test with params file
  #-----------------------------------------------
  
  if [ ${parallel[i]} != "restart" ]; then  
      if [ ${parallel[i]} == "serial" ]; then
        ./flusi ${params_files[i]}
      else
        mpirun -host localhost -n 4 ./flusi ${params_files[i]}
      fi
  else # restart testing needs different commands
      echo -e ${Purple} "-----------------------------" ${Color_Off}
      echo -e ${Purple} "Testing restart capability" ${Color_Off}
      echo -e ${Purple} "-----------------------------" ${Color_Off}
      sleep 0.5

      mpirun -host localhost -n 4 ./flusi Testing_Sphere_start.ini
      echo -e ${Purple} "Waiting 1 sec before restarting..." $Color_Off
      sleep 1
      mpirun -host localhost -n 4 ./flusi Testing_Sphere_restart.ini    
  fi
  
  
  #-----------------------------------------------
  # now we have the new data lying around
  #-----------------------------------------------

  
  for p in ${prefixes[@]}
  do  
    for t in ${times[@]}
    do
      # *.h5 file coming out of the code
      file=${p}"_"${t}".h5"
      # will be transformed into this *.key file
      keyfile=${p}"_"${t}".key"
      # which we will compare to this *.ref file
      reffile=${dirs[i]}"/"${p}"_"${t}".ref"      
      
      if [ -f ${file} ]; then
        echo "comparing " $file " with " $keyfile " and " $reffile

        ./flusi --postprocess --keyvalues ${file}
        # we generated the key file and delete the source now
        rm ${file}
        
        if [ -f $reffile ]; then        
            ./flusi --postprocess --compare-keys $keyfile $reffile            
            if [ $? == "0" ]; then
              echo -e ${Green} "Passed!" ${Color_Off}
              passed=$((passed+1))
            elif [ $? == "1" ]; then
              echo -e ${Red} "ERROR! file=" ${file} " Comparison FAILED!" ${Color_Off}
              echo "It takes me a moment to digest this..."
              sleep 1
              error=$((error+1))
            elif [ $? == "2" ]; then
              echo -e ${Yellow} "WARNING! file=" ${file} " Comparison is not perfect." ${Color_Off}
              echo "It takes me a moment to digest this..."
              sleep 1
              warn=$((warn+1))
            fi
            
            checked=$((checked+1))            
        else
            echo -e ${Red} "Reference file not found" ${Color_Off}            
        fi
        echo "-------------------"
      fi
    done    
  done
  

  #-----------------------------------------------
  # remove the data generated now
  #-----------------------------------------------  
  rm -f *.key
  rm -f *.h5
  rm -f drag_data
done





echo "-------------------------------------"
echo "summary flusi fsi unit testing"
echo "-------------------------------------"
echo -e $Green "Passed: "$passed"/"$checked $Color_Off
if [ $error != 0 ]; then
echo -e $Red "Fail: "$error"/"$checked $Color_Off
fi
if [ $warn != 0 ]; then
echo -e $Yellow "Warned: "$warn"/"$checked $Color_Off
fi
if [ $passed == $checked ]; then
echo -e $Green "All tests passed -> This is good news!" $Color_Off
fi
echo "-------------------------------------"
