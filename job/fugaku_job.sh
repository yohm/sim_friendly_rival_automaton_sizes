#!/bin/bash
#PJM --rsc-list "node=1"
#PJM --rsc-list "rscgrp=small"
#PJM --rsc-list "elapse=00:10:00"
#PJM --mpi "max-proc-per-node=4"
#PJM -s

export OMP_NUM_THREADS=1

mpiexec -stdout-proc ./%/1000R/stdout -stderr-proc ./%/1000R/stderr ../main.out ../friendly_rival_10k
