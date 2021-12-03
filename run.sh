#!/bin/bash


cores=8
# cores=16

# Encourages idle threads to spin rather than sleep.
# export OMP_WAIT_POLICY='active'
# Don't let the runtime deliver fewer threads than those we asked for.
# export OMP_DYNAMIC='false'
export OMP_NUM_THREADS="${cores}"
export GOMP_CPU_AFFINITY="0-$((${cores}-1))"


params=(
    # '207000 8.5 0.84 normal random 1.0 14'
    # '207000 4 0.84 normal random 1.0 14'
    # '207000 4 0.84 normal random 0.2 14'
    '207000 4 0.84 normal diagonal 0.2 14'
    # '207000 20 0.84 normal diagonal 0.5 14'

    # '4000 5 2 normal diagonal 0.05 14'
    # '4000 3.1 125 normal diagonal 0.05 14'
    # '38000 2600 0.43 normal diagonal 0.05 14'
    # '38000 5 125 normal diagonal 0.005 14'
    # '5558000 22 3 normal diagonal 0.05 14'
    # '121000 420 13 normal diagonal 0.05 14'
    # '5558000 100 13 normal random 1 14'

    # '161000 4 0.64 normal diagonal 0.005 14'
    # '161000 10 0.02 normal random 1 14'
    # '207000 4 0.84 normal diagonal 0.05 14'
    # '141000 20 0.64 normal random 1 14'
    # '38000 220 0.84 normal diagonal 0.005 14'
    # '526000 70 2.5 normal diagonal 0.5 14'
    # '4690000 4 1.6 normal random 1 14'
    # '1635000 45 0.2 normal diagonal 0.05 14'
    # '5558000 20 0.64 normal diagonal 0.005 14'

    # '4000 105 0.02 normal random 1 14'
    # '4000 105 0.02 normal diagonal 0.5 14'
    # '4000 105 0.02 normal diagonal 0.05 14'
    # '4000 105 0.02 normal diagonal 0.005 14'

    # '4000 3.1 125 normal random 1 14'
    # '4000 3.1 125 normal diagonal 0.5 14'
    # '4000 3.1 125 normal diagonal 0.05 14'
    # '4000 3.1 125 normal diagonal 0.005 14'

)

# "synthetic" distribution placement diagonal_factor seed nr_rows  nr_cols  nr_nzeros density mem_footprint mem_range  avg_nnz_per_row std_nnz_per_row avg_bw    std_bw    avg_sc  std_sc       format_name  time     gflops W_avg J_estimated
# synthetic,  normal,      random,   1,              14,  4690000, 4690000, 16447528, 0.0001, 143.376,      [128-256], 3.50694,        1.60533,        0.494469, 0.283949, 356223, 1.24247e+06, USE_MKL_IE,  1.76871, 0,     250,  442.176
# rows=4690000, cols=4690000, nnz=16451137, avg_nnz_per_row=3.50771, std_nnz_per_row=1.60551, avg_bw=0.494433, std_bw=0.283955, avg_sc=19.1105, std_sc=6310.95, 

for p in "${params[@]}"; do
    echo $p
    ./artificial_matrix.exe $p
    echo
done



