#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <omp.h>

#include "artificial_matrix_generation.h"


int
main(int argc, char **argv)
{
	struct csr_matrix * csr;
	long nr_rows;
	long nr_cols;
	double avg_nnz_per_row, std_nnz_per_row;
	unsigned int seed;
	char * distribution;
	char * placement;
	double bw;
	double skew;
	double avg_num_neighbours;
	double cross_row_similarity;
	char * matrix_name = "";
	long i;

	if (argc < 6)
	{
		printf("wrong number of parameters\n");
		exit(1);
	}

	i = 1;
	nr_rows = atoi(argv[i++]);
	nr_cols = atoi(argv[i++]);
	avg_nnz_per_row = atof(argv[i++]);
	std_nnz_per_row = atof(argv[i++]);
	distribution = argv[i++];
	placement = argv[i++];
	bw = atof(argv[i++]);
	skew = atof(argv[i++]);
	avg_num_neighbours = atof(argv[i++]);
	cross_row_similarity = atof(argv[i++]);
	seed = atoi(argv[i++]);
	if (argc > i)
		matrix_name = argv[i++];
	else
		matrix_name = "unnamed";
	printf("matrix: %s\n", matrix_name);


	clock_t time_s, time_e;
	time_s = clock();
	csr = artificial_matrix_generation(nr_rows, nr_cols, avg_nnz_per_row, std_nnz_per_row, distribution, seed, placement, bw, skew, avg_num_neighbours, cross_row_similarity);
	time_e = clock();
	printf("time generate matrix = %g\n", ((double) (time_e - time_s))  / CLOCKS_PER_SEC);
	printf("\n---\n\n");
	printf("%s\n", matrix_name);
	printf("distribution = %s\n", csr->distribution);
	printf("placement = %s\n", csr->placement);
	printf("seed = %d\n", csr->seed);

	printf("rows = %d\n", csr->nr_rows);
	printf("cols = %d\n", csr->nr_cols);
	printf("nnz = %d\n", csr->nr_nzeros);

	printf("density = %g\n", csr->density);

	printf("mem_footprint = %g\n", csr->mem_footprint);
	printf("mem_range = %s\n", csr->mem_range);

	printf("avg_nnz_per_row = %g\n", csr->avg_nnz_per_row);
	printf("std_nnz_per_row = %g\n", csr->std_nnz_per_row);
	printf("min_nnz_per_row = %g\n", csr->min_nnz_per_row);
	printf("max_nnz_per_row = %g\n", csr->max_nnz_per_row);

	printf("avg_bw = %g\n", csr->avg_bw);
	printf("std_bw = %g\n", csr->std_bw);
	printf("min_bw = %g\n", csr->min_bw);
	printf("max_bw = %g\n", csr->max_bw);

	printf("avg_bw_scaled = %g\n", csr->avg_bw_scaled);
	printf("std_bw_scaled = %g\n", csr->std_bw_scaled);
	printf("min_bw_scaled = %g\n", csr->min_bw_scaled);
	printf("max_bw_scaled = %g\n", csr->max_bw_scaled);

	printf("avg_sc = %g\n", csr->avg_sc);
	printf("std_sc = %g\n", csr->std_sc);
	printf("min_sc = %g\n", csr->min_sc);
	printf("max_sc = %g\n", csr->max_sc);

	printf("avg_sc_scaled = %g\n", csr->avg_sc_scaled);
	printf("std_sc_scaled = %g\n", csr->std_sc_scaled);
	printf("min_sc_scaled = %g\n", csr->min_sc_scaled);
	printf("max_sc_scaled = %g\n", csr->max_sc_scaled);

	printf("skew = %g\n", csr->skew);

	printf("avg_num_neighbours = %g\n", csr->avg_num_neighbours);
	printf("std_num_neighbours = %g\n", csr->std_num_neighbours);
	printf("min_num_neighbours = %g\n", csr->min_num_neighbours);
	printf("max_num_neighbours = %g\n", csr->max_num_neighbours);

	printf("cross_row_similarity = %g\n", csr->cross_row_similarity);
	// printf("avg_cross_row_similarity = %g\n", csr->avg_cross_row_similarity);
	printf("std_cross_row_similarity = %g\n", csr->std_cross_row_similarity);
	printf("min_cross_row_similarity = %g\n", csr->min_cross_row_similarity);
	printf("max_cross_row_similarity = %g\n", csr->max_cross_row_similarity);
	
	printf("avg_ngroups = %g\n", csr->avg_ngroups);
	printf("std_ngroups = %g\n", csr->std_ngroups);
	printf("min_ngroups = %g\n", csr->min_ngroups);
	printf("max_ngroups = %g\n", csr->max_ngroups);

	printf("avg_ngroups_size = %g\n", csr->avg_ngroups_size);
	printf("std_ngroups_size = %g\n", csr->std_ngroups_size);
	printf("min_ngroups_size = %g\n", csr->min_ngroups_size);
	printf("max_ngroups_size = %g\n", csr->max_ngroups_size);

	printf("avg_clustering = %g\n", csr->avg_clustering);
	printf("std_clustering = %g\n", csr->std_clustering);
	printf("min_clustering = %g\n", csr->min_clustering);
	printf("max_clustering = %g\n", csr->max_clustering);

	printf("avg_dis = %g\n", csr->avg_dis);
	printf("std_dis = %g\n", csr->std_dis);
	printf("min_dis = %g\n", csr->min_dis);
	printf("max_dis = %g\n", csr->max_dis);

	printf("\n");

	// printf("target: %ld\t%ld\t%lf\t%lf\t%s\t%s\t%lf\t%lf\t%lf\t%lf\n", nr_rows, nr_cols, avg_nnz_per_row, std_nnz_per_row, distribution, placement, bw, skew, avg_num_neighbours, cross_row_similarity);
	// printf("result: \t%d\t%d\t%lf\t%lf\t%s\t%s\t%lf\t%lf\t%lf\t%lf\n", csr->nr_rows, csr->nr_cols, csr->avg_nnz_per_row, csr->std_nnz_per_row, csr->distribution, csr->placement, csr->avg_bw_scaled, csr->skew, csr->avg_num_neighbours, csr->cross_row_similarity);

	free_csr_matrix(csr);
	return 0;
}


