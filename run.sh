#!/bin/bash


cores=16

# Encourages idle threads to spin rather than sleep.
# export OMP_WAIT_POLICY='active'
# Don't let the runtime deliver fewer threads than those we asked for.
# export OMP_DYNAMIC='false'
export OMP_NUM_THREADS="${cores}"
export GOMP_CPU_AFFINITY="0-$((${cores}-1))"



./artificial_matrix.exe



