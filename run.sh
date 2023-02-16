#!/bin/bash

cores=4

# Encourages idle threads to spin rather than sleep.
# export OMP_WAIT_POLICY='active'
# Don't let the runtime deliver fewer threads than those we asked for.
# export OMP_DYNAMIC='false'
export OMP_NUM_THREADS="${cores}"
export GOMP_CPU_AFFINITY="0-$((${cores}-1))"

params=(
    '65535 65535 5 1.6667 normal random 0.05 0 0.05 0.05 14'
    '1082401 1082401 10 3.3333 normal random 0.3 10000 0.5 0.5 14'
    '3020732 3020732 50 16.6667 normal random 0.3 100 1.9 0.05 14'
)

for p in "${params[@]}"; do
    ./artificial_matrix.exe $p
    echo
done
