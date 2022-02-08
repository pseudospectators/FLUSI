#!/bin/bash 
#PBS -N unittest
#PBS -q l 
#PBS -b 1 
#PBS -l elapstim_req=2:00:00 
#PBS -l cpunum_job=1 
#PBS -l memsz_job=10gb 

cd $PBS_O_WORKDIR 

export mpi_command="mpirun -f ${NQSII_MPINODES} -ppn 40 -np 4"
export mpi_serial="mpirun -f ${NQSII_MPINODES} -ppn 40 -np 1"

./unittest.sh

# Run as
# qsub job_test.sh
