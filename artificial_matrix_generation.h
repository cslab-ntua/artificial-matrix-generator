#ifndef ARTIFICIAL_MATRIX_GENERATION_H
#define ARTIFICIAL_MATRIX_GENERATION_H


#define ValueType  double


struct csr_matrix {
	int *row_ptr;
	int *col_ind;
	ValueType *values;

	unsigned int nr_rows;
	unsigned int nr_cols;
	unsigned int nr_nzeros;

	double density;
	double mem_footprint;
	int precision;
	char mem_range[128];

	double avg_nnz_per_row;
	double std_nnz_per_row;

	int seed;
	char * distribution;
	char * placement;
	double diagonal_factor;

	double avg_bw;
	double std_bw;
	double avg_sc;
	double std_sc;
};


struct csr_matrix * artificial_matrix_generation(long nr_rows, long nr_cols, double avg_nnz_per_row, double std_nnz_per_row, unsigned int seed, char * placement, double d_f);


#endif /* ARTIFICIAL_MATRIX_GENERATION_H */

