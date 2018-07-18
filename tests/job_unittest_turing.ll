#!/bin/bash
#=========== Global directives ===========
#@ job_name = UNITTEST
#@ shell    = /bin/bash
#============= pre-parallel step =============
#@ step_name = sequential_preprocessing
#@ job_type  = serial
#@ class = archive
#@ cpu_limit = 0:10:00
#@ output = $(job_name).$(jobid).pre
#@ error = $(output)
#@ queue
#============= Parallel step =============
#@ step_name  = parallel_step
#@ dependency = (sequential_preprocessing == 0)
# (submit only if previous step completed without error)
#@ job_type = BLUEGENE
#@ wall_clock_limit = 10:00:00
#@ bg_size = 16
#@ output = $(job_name).$(jobid)
#@ error = $(output)
#@ queue
#============= post-parallel step =============
#@ step_name  = sequential_postprocessing
#@ dependency = (parallel_step >= 0) && (sequential_preprocessing == 0)
# (submit even if previous step completed with an error)
#@ job_type   = serial
#@ class = archive
#@ output = $(job_name).$(jobid).post
#@ error = $(output)
#@ queue


#=========================================
#======= USER-SPECIFIC PART       ========
#=========================================
postprocessing=no
PARAMSFILE=
PROGFILE=flusi
# if dir_gaya is empty, we'll copy to startdir
dir_gaya=takeoff/voluntary_takeoff_video/
files_input_startdir=("*")
files_input_gaya=
RANKS=1
# note $CPU = $RANKS * $bg_size
CPU=4
backup=no
COMMANDS[1]="./$PROGFILE $PARAMSFILE"

cp ../flusi .


case $LOADL_STEP_NAME in
  #=========================================
  #======= Sequential preprocessing ========
  #=========================================
  sequential_preprocessing )
  set -ex
  echo "start dir:" $LOADL_STEP_INITDIR
  if [ ! -f  $LOADL_STEP_INITDIR/$PROGFILE ]; then
    echo  $LOADL_STEP_INITDIR/$PROGFILE "not found!"
    exit
  fi
  # in postprocessing mode, dont care about paramsfile
  if [ "$postprocessing" == "no" ]; then
    if [ ! -f $PARAMSFILE ]; then
      echo  $LOADL_STEP_INITDIR/$PARAMSFILE "not found!"
      exit
    fi
  fi

  # copy all input files to the work directory (from startdir)
  for file in "${files_input_startdir[@]}"
  do
    cp -r $LOADL_STEP_INITDIR/$file $tmpdir
  done

  # copy all input files to the work directory (from server)
  if [ ! "$files_input_gaya" == "" ]; then
    for file in "${files_input_gaya[@]}"
    do
      cd $tmpdir
      mfget $dir_gaya/$file
      cd $LOADL_STEP_INITDIR
    done
  fi


  # symbolic link to workdir
  rm -f workdir
  ln -s $tmpdir workdir
  echo "working dir: " $tmpdir
  cd $tmpdir

  # ---------------------------------------------------------------------------------------------
  # backup resuming goes here
  if [ "$backup" == "yes" ]; then
    if [ "$dir_gaya" == "" ]; then
      echo "no directory on gaya, trying to resume from startdir"
      cp $LOADL_STEP_INITDIR/*.t ./
      cp $LOADL_STEP_INITDIR/runtime* ./
    else
      echo "fetchin input files from ergon..."
      mfget $dir_gaya/*.t
      mfget $dir_gaya/runti*
    fi
  fi
  # ---------------------------------------------------------------------------------------------

  echo "files in starting dir:"
  ls -l
  ;;
  #=========================================
  #============= Parallel step =============
  #=========================================
  parallel_step )
  set -x
  cd $tmpdir
  pwd
  ls -l

  export mpi_command="runjob --np $CPU --ranks-per-node $RANKS --envs BGLOCKLESSMPIO_F_TYPE="0x47504653" : "
  export mpi_serial="runjob --np 1 --ranks-per-node $RANKS --envs BGLOCKLESSMPIO_F_TYPE="0x47504653" : "

  ./unittest.sh

  ;;
  #=========================================
  #======= Sequential postprocessing =======
  #=========================================
  sequential_postprocessing )
  set -x
  # remove link to workdir
  rm workdir
  # goto tempdir and mv all files to gaya
  cd $tmpdir
  echo "files in work dir"
  ls -l

  for file in *.log
  do
    cp $file $LOADL_STEP_INITDIR/
  done
  ;;
esac
