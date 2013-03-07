# Reset
Color_Off='\e[0m'       # Text Reset

# Pretty colors
# Regular Colors
Black='\e[0;30m'        # Black
Red='\e[0;31m'          # Red
Green='\e[0;32m'        # Green
Yellow='\e[0;33m'       # Yellow
Blue='\e[0;34m'         # Blue
Purple='\e[0;35m'       # Purple
Cyan='\e[0;36m'         # Cyan
White='\e[0;37m'        # White

# Bold
BBlack='\e[1;30m'       # Black
BRed='\e[1;31m'         # Red
BGreen='\e[1;32m'       # Green
BYellow='\e[1;33m'      # Yellow
BBlue='\e[1;34m'        # Blue
BPurple='\e[1;35m'      # Purple
BCyan='\e[1;36m'        # Cyan
BWhite='\e[1;37m'       # White

nproc="2"
params="PARAMSs.ring"

echo -e ${BBlue} "deleting main..." ${Color_Off}
rm main

echo -e ${BBlue} "make clean..." ${Color_Off}
make clean

echo -e ${BBlue} "make sugiton..." ${Color_Off}
make sugiton

if [ -f main ]
then
cp main Test/
cd Test/
rm *mpiio*
echo "-----------"
echo -e ${BCyan} "mpirun -n" ${nproc} "./main <" ${params} ${Color_Off}
echo "-----------"
mpirun -n ${nproc} ./main < ${params}

else
echo -e ${BRed} "MAKE FAILED!" ${Color_Off}
fi
