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

	double skew;

	int seed;
	char * distribution;
	char * placement;

	double avg_nnz_per_row;
	double std_nnz_per_row;
	double min_nnz_per_row;
	double max_nnz_per_row;

	double avg_sc;
	double std_sc;
	double min_sc;
	double max_sc;
	
	double avg_sc_scaled;
	double std_sc_scaled;
	double min_sc_scaled;
	double max_sc_scaled;

	double avg_bw;
	double std_bw;
	double min_bw;
	double max_bw;

	double avg_bw_scaled;
	double std_bw_scaled;
	double min_bw_scaled;
	double max_bw_scaled;

	double avg_num_neighbours;
	double std_num_neighbours;
	double min_num_neighbours;
	double max_num_neighbours;

	double cross_row_similarity;
	// double avg_cross_row_similarity;
	double std_cross_row_similarity;
	double min_cross_row_similarity;
	double max_cross_row_similarity;

	double avg_clustering;
	double std_clustering;
	double min_clustering;
	double max_clustering;

	double avg_dis;
	double std_dis;
	double min_dis;
	double max_dis;

	double avg_ngroups;
	double std_ngroups;
	double min_ngroups;
	double max_ngroups;

	double avg_ngroups_size;
	double std_ngroups_size;
	double min_ngroups_size;
	double max_ngroups_size;

};


int free_csr_matrix(struct csr_matrix * csr);
struct csr_matrix * artificial_matrix_generation(long nr_rows, long nr_cols, double avg_nnz_per_row, double std_nnz_per_row, char * distribution, unsigned int seed, char * placement, double bw_scaled, double skew, double avg_num_neighbours, double cross_row_similarity);
void csr_matrix_print(struct csr_matrix * csr);
void csr_matrix_write_mtx(struct csr_matrix * csr, char * file_out);


#endif /* ARTIFICIAL_MATRIX_GENERATION_H */

