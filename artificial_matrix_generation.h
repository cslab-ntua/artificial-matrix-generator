#ifndef ARTIFICIAL_MATRIX_GENERATION_H
#define ARTIFICIAL_MATRIX_GENERATION_H


#ifndef ValueType
	#define ValueType  double
#endif


struct csr_matrix {
	int * row_ptr;
	int * col_ind;
	ValueType * values;

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
	double bandwidth_scaled;

	double avg_bw;
	double std_bw;
	double avg_sc;
	double std_sc;
};


int free_csr_matrix(struct csr_matrix * csr);
struct csr_matrix * artificial_matrix_generation(long nr_rows, long nr_cols, double avg_nnz_per_row, double std_nnz_per_row, char * distribution, unsigned int seed, char * placement, double bw_scaled);
void csr_matrix_print(struct csr_matrix * csr);
void csr_matrix_write_mtx(struct csr_matrix * csr, char * file_out);


#endif /* ARTIFICIAL_MATRIX_GENERATION_H */

