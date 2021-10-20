#!/bin/bash


cores=16

# Encourages idle threads to spin rather than sleep.
# export OMP_WAIT_POLICY='active'
# Don't let the runtime deliver fewer threads than those we asked for.
# export OMP_DYNAMIC='false'
export OMP_NUM_THREADS="${cores}"
export GOMP_CPU_AFFINITY="0-$((${cores}-1))"


params=(

    # '4000 5 2 normal diagonal 0.05 14'
    # '4000 3.1 125 normal diagonal 0.05 14'
    # '38000 2600 0.43 normal diagonal 0.05 14'
    # '38000 5 125 normal diagonal 0.005 14'
    '5558000 22 3 normal diagonal 0.05 14'


    # '4000 105 0.02 normal random 1 14'
    # '4000 105 0.02 normal diagonal 0.5 14'
    # '4000 105 0.02 normal diagonal 0.05 14'
    # '4000 105 0.02 normal diagonal 0.005 14'

    # '4000 3.1 125 normal random 1 14'
    # '4000 3.1 125 normal diagonal 0.5 14'
    # '4000 3.1 125 normal diagonal 0.05 14'
    # '4000 3.1 125 normal diagonal 0.005 14'

)


for p in "${params[@]}"; do
    ./artificial_matrix.exe $p
done



