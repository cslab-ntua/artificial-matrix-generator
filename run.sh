#!/bin/bash


cores=24
# cores=16

# Encourages idle threads to spin rather than sleep.
# export OMP_WAIT_POLICY='active'
# Don't let the runtime deliver fewer threads than those we asked for.
# export OMP_DYNAMIC='false'
export OMP_NUM_THREADS="${cores}"
export GOMP_CPU_AFFINITY="0-$((${cores}-1))"


distribution='normal'
# distribution='gamma'


params=(

    # "170998 170998 5.607878455 4.39216211 ${distribution} random 0.2972525308 61.94715601 0.803365 0.633019 14 scircuit"
    # "206500 206500 6.166532688 4.435865332 ${distribution} random 0.001907956522 6.135290158 0.176686 0.330509 14 mac_econ_fwd500"
    # "21200 21200 70.22490566 6.326999832 ${distribution} random 0.0662003204 0.1391969736 1.916 0.963020 14 raefsky3"
    # "38744 38744 45.72893867 38.39531373 ${distribution} random 0.02989415739 1.755366813 1.262626 0.853730 14 bbmat"
    # "49152 49152 39 0 ${distribution} random 0.2446905772 0 1.4415068 0.810948 14 conf5_4-8x-15"
    # "525825 525825 3.994152047 0.07632277052 ${distribution} random 0.001342169398 0.001464128789 0.498297 0.998906 14 mc2depi"
    # "46835 46835 50.68860895 27.78059606 ${distribution} random 0.1877712645 1.86060326 1.719725 0.866408 14 rma10"
    # "121192 121192 21.65432537 13.79266245 ${distribution} random 0.6230549795 2.740592173 1.095827 0.633400 14 cop20k_A"
    # "1000005 1000005 3.105520472 25.34520973 ${distribution} random 0.1524629001 1512.433913 0.315681 0.936095 14 webbase-1M"
    # "62451 62451 64.16843605 14.05626099 ${distribution} random 0.008604097649 0.215550897 1.61575 0.914729 14 cant"
    # "36417 36417 119.305956 31.86038422 ${distribution} random 0.1299377034 0.7098894878 1.837758 0.931726 14 pdb1HYS"
    # "42138 42138 104.73798 102.4431672 ${distribution} random 0.6066963064 0.9954557077 1.923459 0.992154 14 TSOPF_RS_b300_c3"
    # "68121 68121 78.94424627 1061.43997 ${distribution} random 0.04535742434 861.9001253 0.895986 0.993706 14 Chebyshev4"
    # "83334 83334 72.125183 19.08019415 ${distribution} random 0.06981133797 0.1230474105 1.71335 0.882633 14 consph"
    # "140874 140874 55.46377614 11.07481064 ${distribution} random 0.04587554507 0.8390381452 1.712385 0.873444 14 shipsec1"
    # "161070 161070 50.81725958 19.6982847 ${distribution} random 0.03958570835 0.8104085258 1.279542 0.868860 14 PR02R"
    # "66463 66463 155.7681567 350.7443175 ${distribution} random 0.5902849043 425.2424452 1.924902 0.961470 14 mip1"
    # "4284 1096894 2633.994398 4209.259315 ${distribution} random 0.9555418052 20.32958219 1.576267 0.269166 14 rail4284"
    # "217918 217918 53.38899953 4.743895102 ${distribution} random 0.05932070192 2.371481046 1.878352 0.949809 14 pwtk"
    # "63838 63838 221.6369247 95.87570167 ${distribution} random 0.865768503 14.44417747 1.731501 0.866870 14 crankseg_2"
    # "185639 185639 80.86266894 126.9718576 ${distribution} random 0.1935377211 7.186719641 1.269782 0.963174 14 Si41Ge41H72"
    # "38120 38120 424.2174449 484.237499 ${distribution} random 0.4832654419 1.317207865 1.979904 0.927411 14 TSOPF_RS_b2383"
    # "1382908 1382908 12.23295621 37.23001003 ${distribution} random 0.02150049625 632.7797558 1.507183 0.768973 14 in-2004"
    # "268096 268096 68.96214789 105.3875532 ${distribution} random 0.1722725043 9.179497325 1.184016 0.969575 14 Ga41As41H72"
    # "862664 862664 22.29737186 29.33341123 ${distribution} random 0.248901028 312.2656191 1.330828 0.744158 14 eu-2005"
    # "1634989 1634989 12.08147455 31.07497821 ${distribution} random 0.3400604749 410.3736266 0.067098 0.069594 14 wikipedia-20051105"
    # "4690002 4690002 4.33182182 1.106157847 ${distribution} random 0.0009716598294 288.0238916 0.553744 0.648584 14 rajat31"
    # "952203 952203 48.85772782 11.94657153 ${distribution} random 0.2042067138 0.5760045224 1.79674 0.906047 14 ldoor"
    # "5558326 5558326 10.709032 1356.616274 ${distribution} random 0.5025163836 120504.8496 1.065325 0.942768 14 circuit5M"
    # "986703 986703 72.63211422 15.81042955 ${distribution} random 0.01817301802 0.1152091726 1.769845 0.915136 14 bone010"
    # "5154859 5154859 19.24389222 5.736719369 ${distribution} random 0.2119567644 1.442333363 0.197548 0.794106 14 cage15"

    # "986703 986703 72.63211422 0 ${distribution} random 0.0922  3 1.5577 0.6427 14 synthetic"
    # "986703 986703 72.63211422 1 ${distribution} random 0.0922  3 1.5577 0.6427 14 synthetic"
    # "986703 986703 72.63211422 2 ${distribution} random 0.0922  3 1.5577 0.6427 14 synthetic"
    # "986703 986703 72.63211422 4 ${distribution} random 0.0922  3 1.5577 0.6427 14 synthetic"
    # "986703 986703 72.63211422 8 ${distribution} random 0.0922  3 1.5577 0.6427 14 synthetic"
    # "986703 986703 72.63211422 16 ${distribution} random 0.0922 3 1.5577 0.6427 14 synthetic"
    # "986703 986703 72.63211422 32 ${distribution} random 0.0922 3 1.5577 0.6427 14 synthetic"
    # "986703 986703 72.63211422 64 ${distribution} random 0.0922 3 1.5577 0.6427 14 synthetic"

    # "819190 819190 50 1.6667 ${distribution} random 0.05 0 0.05 0.05 14 synthetic"
    # "819190 819190 50 1.6667 ${distribution} random 0.05 0 0.05 0.5 14 synthetic"
    # "819190 819190 50 1.6667 ${distribution} random 0.05 0 0.05 0.95 14 synthetic"
    # "819190 819190 50 1.6667 ${distribution} random 0.05 0 0.5 0.05 14 synthetic"
    # "819190 819190 50 1.6667 ${distribution} random 0.05 0 0.5 0.5 14 synthetic"
    # "819190 819190 50 1.6667 ${distribution} random 0.05 0 0.5 0.95 14 synthetic"
    # "819190 819190 50 1.6667 ${distribution} random 0.05 0 0.95 0.05 14 synthetic"
    # "819190 819190 50 1.6667 ${distribution} random 0.05 0 0.95 0.5 14 synthetic"
    # "819190 819190 50 1.6667 ${distribution} random 0.05 0 0.95 0.95 14 synthetic"
    # "819190 819190 50 1.6667 ${distribution} random 0.05 100 0.05 0.05 14 synthetic"
    # "819190 819190 50 1.6667 ${distribution} random 0.05 100 0.05 0.5 14 synthetic"
    # "819190 819190 50 1.6667 ${distribution} random 0.05 100 0.05 0.95 14 synthetic"
    # "819190 819190 50 1.6667 ${distribution} random 0.05 100 0.5 0.05 14 synthetic"
    # "819190 819190 50 1.6667 ${distribution} random 0.05 100 0.5 0.5 14 synthetic"
    # "819190 819190 50 1.6667 ${distribution} random 0.05 100 0.5 0.95 14 synthetic"
    # "819190 819190 50 1.6667 ${distribution} random 0.05 100 0.95 0.05 14 synthetic"
    # "819190 819190 50 1.6667 ${distribution} random 0.05 100 0.95 0.5 14 synthetic"
    # "819190 819190 50 1.6667 ${distribution} random 0.05 100 0.95 0.95 14 synthetic"
    # "819190 819190 50 1.6667 ${distribution} random 0.05 1000 0.05 0.05 14 synthetic"
    # "819190 819190 50 1.6667 ${distribution} random 0.05 1000 0.05 0.5 14 synthetic"
    # "819190 819190 50 1.6667 ${distribution} random 0.05 1000 0.05 0.95 14 synthetic"
    # "819190 819190 50 1.6667 ${distribution} random 0.05 1000 0.5 0.05 14 synthetic"
    # "819190 819190 50 1.6667 ${distribution} random 0.05 1000 0.5 0.5 14 synthetic"
    # "819190 819190 50 1.6667 ${distribution} random 0.05 1000 0.5 0.95 14 synthetic"
    # "819190 819190 50 1.6667 ${distribution} random 0.05 1000 0.95 0.05 14 synthetic"
    # "819190 819190 50 1.6667 ${distribution} random 0.05 1000 0.95 0.5 14 synthetic"
    # "819190 819190 50 1.6667 ${distribution} random 0.05 1000 0.95 0.95 14 synthetic"
    # "819190 819190 50 1.6667 ${distribution} random 0.05 10000 0.05 0.05 14 synthetic"
    # "819190 819190 50 1.6667 ${distribution} random 0.05 10000 0.05 0.5 14 synthetic"
    # "819190 819190 50 1.6667 ${distribution} random 0.05 10000 0.05 0.95 14 synthetic"
    # "819190 819190 50 1.6667 ${distribution} random 0.05 10000 0.5 0.05 14 synthetic"
    # "819190 819190 50 1.6667 ${distribution} random 0.05 10000 0.5 0.5 14 synthetic"
    # "819190 819190 50 1.6667 ${distribution} random 0.05 10000 0.5 0.95 14 synthetic"
    # "819190 819190 50 1.6667 ${distribution} random 0.05 10000 0.95 0.05 14 synthetic"
    # "819190 819190 50 1.6667 ${distribution} random 0.05 10000 0.95 0.5 14 synthetic"
    # "819190 819190 50 1.6667 ${distribution} random 0.05 10000 0.95 0.95 14 synthetic"

    # "38120 38120 424.2174449 484.237499 ${distribution} random 0.4832654419 1.317207865 1.979904 14 TSOPF_RS_b2383"
    # "38120 38120 424.2174449 484.237499 ${distribution} random 0.4832654419 1.317207865 1 14 TSOPF_RS_b2383"
    # "38120 38120 424.2174449 484.237499 ${distribution} random 0.4832654419 1.317207865 0 14 TSOPF_RS_b2383"

    # '28508159 28508159 5 1.6667 normal random 0.05 0 0.05 0.5 14'

    # '303884 303884 500 166.6667 normal random 0.05 1000 0.05 0.05 14'

    # '65535 65535 5 1.6667 normal 14 random 0.05 0 0.05 0.05'

    '65535 65535 5 1.6667 normal random 0.05 0 0.05 0.05 14'
    '1082401 1082401 10 3.3333 normal random 0.3 10000 0.5 0.5 14'
    '3020732 3020732 50 16.6667 normal random 0.3 100 1.9 0.05 14'

    # '256 1 1 0 normal random 0.5 0 0 0 14 test'

)

# "synthetic" distribution placement diagonal_factor seed nr_rows  nr_cols  nr_nzeros density mem_footprint mem_range  avg_nnz_per_row std_nnz_per_row avg_bw    std_bw    avg_sc  std_sc       format_name  time     gflops W_avg J_estimated
# synthetic,  normal,      random,   1,              14,  4690000, 4690000, 16447528, 0.0001, 143.376,      [128-256], 3.50694,        1.60533,        0.494469, 0.283949, 356223, 1.24247e+06, USE_MKL_IE,  1.76871, 0,     250,  442.176
# rows=4690000, cols=4690000, nnz=16451137, avg_nnz_per_row=3.50771, std_nnz_per_row=1.60551, avg_bw=0.494433, std_bw=0.283955, avg_sc=19.1105, std_sc=6310.95, 

for p in "${params[@]}"; do
    ./artificial_matrix.exe $p
    echo
done



