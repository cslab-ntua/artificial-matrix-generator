# artificial-matrix-generator
This project generates artificial matrices, given specific structural features (parameters), than can be used for several sparse-matrix benchmarks. The features that have to be provided to the generator are the following : 

* number of rows
* number of columns
* average nonzeros per row
* standard deviation of nonzeros per row
* distribution of nonzeros per row (currently "normal" and "gamma" distributions are supported)
* seed for the random number generator (need to reproduce results on the same artificial matrices across different platforms)
* placement of nonzeros within row (currently "random", "diagonal" and "simple" are supported)
* bandwidth of matrix, which confines the column range that a nonzero can be placed at
* coefficient of skewness, which indicates the imbalance of row length
* average number of neighboring nonzeros within row (neighbors : maximum distance 1 for column values)
* similarity between adjacent rows (based on column position of nonzeros)

The output is a struct with the CSR-representation (row-pointer, column-index and values) of the generated sparse matrix and the above mentioned structural features.

An example of the usage of the generator can be found at [`artificial_matrix.c`](./artificial_matrix.c) (called via bash script [`run.sh`](./run.sh))
