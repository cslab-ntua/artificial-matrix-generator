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
	printf("time generate matrix = %g\n", ((double) (time_s - time_e))  / CLOCKS_PER_SEC);

	printf("%s, ", matrix_name);
	printf("distribution=%s, ", csr->distribution);
	printf("placement=%s, ", csr->placement);
	printf("seed=%d, ", csr->seed);
	printf("rows=%d, ", csr->nr_rows);
	printf("cols=%d, ", csr->nr_cols);
	printf("nnz=%d, ", csr->nr_nzeros);
	printf("density=%g, ", csr->density);
	printf("mem_footprint=%g, ", csr->mem_footprint);
	printf("mem_range=%s, ", csr->mem_range);
	printf("avg_nnz_per_row=%g, ", csr->avg_nnz_per_row);
	printf("std_nnz_per_row=%g, ", csr->std_nnz_per_row);
	printf("avg_bw=%g, ", csr->avg_bw);
	printf("std_bw=%g, ", csr->std_bw);
	printf("avg_bw_scaled=%g, ", csr->avg_bw_scaled);
	printf("std_bw_scaled=%g, ", csr->std_bw_scaled);
	printf("avg_sc=%g, ", csr->avg_sc);
	printf("std_sc=%g, ", csr->std_sc);
	printf("avg_sc_scaled=%g, ", csr->avg_sc_scaled);
	printf("std_sc_scaled=%g, ", csr->std_sc_scaled);
	printf("min_nnz_per_row=%g, ", csr->min_nnz_per_row);
	printf("max_nnz_per_row=%g, ", csr->max_nnz_per_row);
	printf("skew=%g, ", csr->skew);
	printf("avg_num_neighbours=%g, ", csr->avg_num_neighbours);
	printf("cross_row_similarity=%g, ", csr->cross_row_similarity);
	printf("\n");

	printf("target: %ld\t%ld\t%lf\t%lf\t%s\t%s\t%lf\t%lf\t%lf\t%lf\n", nr_rows, nr_cols, avg_nnz_per_row, std_nnz_per_row, distribution, placement, bw, skew, avg_num_neighbours, cross_row_similarity);
	printf("result: \t%d\t%d\t%lf\t%lf\t%s\t%s\t%lf\t%lf\t%lf\t%lf\n", csr->nr_rows, csr->nr_cols, csr->avg_nnz_per_row, csr->std_nnz_per_row, csr->distribution, csr->placement, csr->avg_bw_scaled, csr->skew, csr->avg_num_neighbours, csr->cross_row_similarity);

	free_csr_matrix(csr);
	return 0;
}


