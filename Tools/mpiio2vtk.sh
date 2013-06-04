#!/bin/bash
#-------------------------------------------------------------------------------
# convert files from mpi format to a PARAVIEW readable format
#-------------------------------------------------------------------------------

# usage: 
# mpiio2vtk.script <NX> <NY> <NZ>
# where NX, NY, and NZ are the number of elements in each dimension

USAGE="mpiio2vtk.script <NX> <NY> <NZ>"

if [ "$1" == "" ]; then
    echo "need to specify NX"
    echo $USAGE
    exit 1
fi
if [ "$2" == "" ]; then
    echo "need to specify NY"
    echo $USAGE
    exit 1
fi
if [ "$3" == "" ]; then
    echo "need to specify NZ"
    echo $USAGE
    exit 1
fi

NX=$1
NY=$2
NZ=$3

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

echo -e $Green "****************************" $Color_Off
echo -e $Green "**      mpiio2vtk         **" $Color_Off
echo -e $Green "****************************" $Color_Off


# make directory
echo "this script will use the subdirectory ./vtk/ asput directory."
echo -e $Cyan "any key to continue" $Color_Off
read dummy

ask_to_skip="yes"
if [ -d vtk/ ]
then
    echo "subdirectory already exists"
    echo -e $Cyan "do you want to clean it? (y/[n])" $Color_Off
    read ask
    if [ "$ask" == "y" ]
    then
	echo -e ${Red} "deleting" ${Color_Off}
	rm -r vtk
	mkdir vtk
	ask_to_skip="no"
    fi  
else
    echo "creating subdirectory..."
    mkdir vtk
    ask_to_skip="no"
fi


# ask what to do with existing *.vtk files
if [ "${ask_to_skip}" == "yes" ]; then
    echo -e ${Cyan} "If a vtk file already exists, do you want to skip it? ([y],n)" ${Color_Off}
    read skip
    if [ "${skip}" == "" ]; then
	skip="y"
    fi
else
    skip="y"
fi

# ask what modus to use
echo "You can either save all variables to ONE *.vtk file"
echo "or you can save each scalar to one file and all vector components to one file"
echo -e ${Cyan} "Do you want to use the several vtk-files mode? ([y],n)" ${Color_Off}
read modus
if [ "${modus}" == "" ]; then
    modus="y"
fi

# ask what file ending to expect
echo -e ${Cyan} "What type of source files to you want to convert, *.binary or *.mpiio? (b,[m])" ${Color_Off}
read ending
if [ "${ending}" == "b" ]; then
    ending="binary"
else
    ending="mpiio"
fi

# find prefixes
# look through all the files whose names end with *.mpiio and put
# them in the items array, as well as in the list, where file names
# are separated with colons. N is the number of items in the list.
N=0
lastp=""
for F in `ls *.${ending}`
do
    # as is the file name with everything after 
    p=$(echo ${F}  | sed 's/_[^_]*$//')_ 
    p=${p%%_}
    if [ "$p" != "$lastp" ] ; then
	lastp=$p
	items[$N]=$p
	N=$((N+1))
    fi
done

echo -e "prefixes  : " ${Blue} ${items[@]} ${Color_Off}


# Begin conversion

if [ "${modus}" == "n" ]
then 
    echo -e ${Cyan} "-----------------------" ${Color_Off}
    echo -e ${Cyan} "ONE vtk file modus" ${Color_Off}
    echo -e ${Cyan} "-----------------------" ${Color_Off}
    lead_prefix=${items[0]}"_"
    
    for file in ${lead_prefix}*.${ending} 
# note: I suppose that you saved the mask
    do
	base=${file%%.${ending}}   # remove trailing .mpiio
	base=${base##${lead_prefix}}    # remove leading ux_

    # note program also works if one or more files are not present
	if [ ${skip}=="y" ]
	then      
	    if [ -f ./vtk/para_${base}.vtk ]
	    then
		echo "file ./vtk/para_"${base}".vtk exists. skipping..."
	    else
 		./convert_mpiio2vtk_ALL ./vtk/para_${base}.vtk ${NX} ${NY} ${NZ} "ux_"${base}"."${ending} "uy_"${base}"."${ending} "uz_"${base}"."${ending} "vorx_"${base}"."${ending} "vory_"${base}"."${ending} "vorz_"${base}"."${ending} "p_"${base}"."${ending} "mask_"${base}"."${ending} "usx_"${base}"."${ending} "usy_"${base}"."${ending} "usz_"${base}"."${ending}
	    fi      
	else
 	    ./convert_mpiio2vtk_ALL ./vtk/para_${base}.vtk ${NX} ${NY} ${NZ} "ux_"${base}"."${ending} "uy_"${base}"."${ending} "uz_"${base}"."${ending} "vorx_"${base}"."${ending} "vory_"${base}"."${ending} "vorz_"${base}"."${ending} "p_"${base}"."${ending} "mask_"${base}"."${ending} "usx_"${base}"."${ending} "usy_"${base}"."${ending} "usz_"${base}"."${ending}
	fi    
    done

else 


    echo -e ${Cyan} "-----------------------" ${Color_Off}
    echo -e ${Cyan} "SEVERAL vtk files modus" ${Color_Off}
    echo -e ${Cyan} "-----------------------" ${Color_Off}
    
    # indentify vectors and scalars from the prefixes

    # look through all prefixed in array items[] if a prefix ends with
    # "x", it's assumed to be the x-component of a vector.  the next 2
    # indices then will contain "y" and "z" component. we remove them
    # from the list items[] the prefix ending with "x" is then
    # converted to its root (remove the trailing "x") and added to the
    # list of vectors[] otherwise, we add the prefix to the list of
    # scalars[]
    N2=0
    N3=0
    for (( i=0; i<N; i++ ))
    do
      # the prefix
	p=${items[i]}
	if [ "${p:${#p}-1:${#p}}" == "x" ]; then 
        # is the last char an "x"? yes -> delete following entries
	    unset items[i+1] 
# delete next two entrys (with y and z ending, hopefully)
	    unset items[i+2]    
	# the trailing "x" indicates a vector
	    vectors[N2]=${p%%x} # collect entries for vector fields
	    N2=$((N2+1))
	else
        # no? it's a scalar.
            if [ "$p" != "" ]; then  # note empty values are not scalars (they are uy, uz but unset because of ux)
		scalars[N3]=${p}
		N3=$((N3+1))
	    fi
	fi
    done
    
    # print summary
    echo -e "found scalars: " ${Cyan} ${scalars[@]} ${Color_Off}
    echo -e "found vectors: " ${Cyan} ${vectors[@]} ${Color_Off}
    echo -e $Cyan "any key to continue" $Color_Off
    read dummy
    
    # --------------------------
    # SCALARS
    # --------------------------
    for (( i=0; i<N3; i++ ))
    do
      # the prefix
	p=${scalars[i]}
      # find the files that start with the prefix
	FLIST=$( ls ${p}*.${ending} )

      # loop over the files starting with the prefix and process
	for F in ${FLIST}
	do
	    vtk=${F%%${ending}}vtk
	  # check if target file exists in vtk/
	    if [ -f vtk/${vtk} ]; then
	      # if you shouldn't skip it, overwrite it.
		if [ "${skip}" != "y" ]; then
		    echo -e $Green "file vtk/"${vtk} "exists, overwriting." $Color_Off
		    convert_mpiio2vtk ${NX} ${NY} ${NZ} -scalar ${F}
		# currently, convert_mpiio2vtk creates the vtk file in the same directory. so move it.
		    mv ${vtk} vtk/
		else
		# file exists and no overwriting.
		    echo -e $Green "file vtk/"${vtk} "exists, skipping." $Color_Off
		fi
	    else
	    # file does not exist, go for it.
		./convert_mpiio2vtk ${NX} ${NY} ${NZ} -scalar ${F}
	    # currently, convert_mpiio2vtk creates the vtk file in the same directory. so move it.
		mv ${vtk} vtk/  
	    fi
	done
    done
    
    # VECTORS
    for (( i=0; i<N2; i++ ))
    do
      # the prefix
	p=${vectors[i]}
      # find the files that start with the prefix + x-component
	FLIST=$( ls ${p}x*.${ending} )

      # loop over the files starting with the prefix and process with raw2vdf
	for F in ${FLIST}
	do
	    vtk=${p}${F##${p}x}  # F=vorx_00010.mpiio vtk=vor_00010.mpiio
	    vtk=${vtk%%${ending}}vtk # vtk=vor_00010.vtk
	    
	  # check if target file exists in vtk/
	    if [ -f vtk/${vtk} ]; then
	      # if you shouldn't skip it, overwrite it.
		if [ "${skip}" != "y" ]; then
		    echo -e $Green "file vtk/"${vtk} "exists, overwriting." $Color_Off
		    ./convert_mpiio2vtk ${NX} ${NY} ${NZ} -vector ${F} ${p}y${F##${p}x} ${p}z${F##${p}x}
		# currently, convert_mpiio2vtk creates the vtk file in the same directory. so move it.
		    mv ${vtk} vtk/
		else
		# file exists and no overwriting.
		    echo -e $Green "file vtk/"${vtk} "exists, skipping." $Color_Off
		fi
	    else
	    # file does not exist, go for it.
		./convert_mpiio2vtk ${NX} ${NY} ${NZ} -vector ${F} ${p}y${F##${p}x} ${p}z${F##${p}x}
	    # currently, convert_mpiio2vtk creates the vtk file in the same directory. so move it.
		mv ${vtk} vtk/
	    fi	    
	done
    done
    
    
    
fi
