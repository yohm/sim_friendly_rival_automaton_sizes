#!/bin/bash
#PJM --rsc-list "node=800"
#PJM --rsc-list "rscgrp=large"
#PJM --rsc-list "elapse=04:00:00"
#PJM --mpi "max-proc-per-node=4"
#PJM -s

export OMP_NUM_THREADS=12

mpiexec -stdout-proc ./%/1000R/stdout -stderr-proc ./%/1000R/stderr ../../main.out dis.passed 9
