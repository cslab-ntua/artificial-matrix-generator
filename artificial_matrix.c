#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

#include "artificial_matrix_generation.h"
#include "matrix_util.h"


#include "debug.h"
#include "time_it.h"
#include "plot/plot.h"


void
plot_csr(struct csr_matrix * csr, char * matrix_name)
{
	long num_pixels = 1024;
	long num_pixels_x = num_pixels, num_pixels_y = num_pixels;
	long buf_n = strlen(matrix_name) + 1 + 1000;
	char buf[buf_n], buf_title[buf_n];

	double get_row(void * A, int pos)
	{
		struct csr_matrix * csr = A;
		long i = binary_search(csr->row_ptr, 0, csr->nr_rows, pos, NULL, NULL);
		if (csr->row_ptr[i] > pos)
			i--;
		return (double) i;
	}
	double get_col(void * A, int pos)
	{
		struct csr_matrix * csr = A;
		return (double) csr->col_ind[pos];
	}
	double get_degree(void * A, int i)
	{
		struct csr_matrix * csr = A;
		return (double) csr->row_ptr[i+1] -  csr->row_ptr[i];
	}

	snprintf(buf, buf_n, "figures/%s-artificial.png", matrix_name);
	snprintf(buf_title, buf_n, "%s-artificial", matrix_name);
	figure_simple_plot(buf, num_pixels_x, num_pixels_y, (csr, csr, NULL, csr->nr_nzeros, 0, get_col, get_row),
		figure_axes_flip_y(_fig);
		figure_enable_legend(_fig);
		figure_set_title(_fig, buf_title);
		figure_set_bounds_x(_fig, 0, csr->nr_cols);
		figure_set_bounds_y(_fig, 0, csr->nr_rows);
	);

	// Plot degree histogram.
	#if 1
		long num_bins = 10000;
		snprintf(buf, buf_n, "figures/%s-artificial_degree_distribution.png", matrix_name);
		snprintf(buf_title, buf_n, "%s-artificial: degree distribution (percentages)", matrix_name);
		figure_simple_plot(buf, num_pixels_x, num_pixels_y, (NULL, csr, NULL, csr->nr_rows, 0, NULL, get_degree),
			figure_enable_legend(_fig);
			figure_set_title(_fig, buf_title);
			figure_series_type_histogram(_s, num_bins, 1);
			figure_series_type_barplot(_s);
		);
	#endif

	// Plot neighbour distances frequencies.
	#if 0
		double * neigh_distances_freq_perc = csr_neighbours_distances_frequencies(csr->row_ptr, csr->col_ind, csr->nr_rows, csr->nr_cols, csr->nr_nzeros);
		snprintf(buf_title, buf_n, "%s-artificial: neighbour distances frequencies (percentages)", matrix_name);
		long pos;
		for (pos=csr->nr_cols-1;pos>0;pos--)
			if (neigh_distances_freq_perc[pos] != 0)
				break;
		snprintf(buf, buf_n, "figures/%s-artificial_neigh_dist_freq.png", matrix_name);
		figure_simple_plot(buf, num_pixels_x, num_pixels_y, (NULL, neigh_distances_freq_perc, NULL, pos+1, 0),
			figure_enable_legend(_fig);
			figure_set_title(_fig, buf_title);
			figure_series_type_barplot(_s);
		);
		snprintf(buf, buf_n, "figures/%s-artificial_neigh_dist_freq_x_100.png", matrix_name);
		figure_simple_plot(buf, num_pixels_x, num_pixels_y, (NULL, neigh_distances_freq_perc, NULL, pos+1, 0),
			figure_enable_legend(_fig);
			figure_set_title(_fig, buf_title);
			figure_series_type_barplot(_s);
			figure_set_bounds_x(_fig, 0, 100);
		);
		snprintf(buf, buf_n, "figures/%s-artificial_neigh_dist_freq_x_1000.png", matrix_name);
		figure_simple_plot(buf, num_pixels_x, num_pixels_y, (NULL, neigh_distances_freq_perc, NULL, pos+1, 0),
			figure_enable_legend(_fig);
			figure_set_title(_fig, buf_title);
			figure_series_type_barplot(_s);
			figure_set_bounds_x(_fig, 0, 1000);
		);
		free(neigh_distances_freq_perc);
	#endif
}


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
		error("wrong number of parameters\n");

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
	if (argc >= i)
	{
		matrix_name = argv[i++];
		printf("matrix: %s\n", matrix_name);
	}
	else
		matrix_name = "csr";


	double time;
	time = time_it(1,
		csr = artificial_matrix_generation(nr_rows, nr_cols, avg_nnz_per_row, std_nnz_per_row, distribution, seed, placement, bw, skew, avg_num_neighbours, cross_row_similarity);
	);
	printf("time generate matrix = %g\n", time);

	// long window_size;  // Distance from left and right.
	// window_size = 64 / sizeof(double) / 2;
	// window_size = 8;
	// window_size = 4;
	// window_size = 1;
	// double avg_num_neigh = csr_avg_row_neighbours(csr->row_ptr, csr->col_ind, csr->nr_rows, csr->nr_cols, csr->nr_nzeros, window_size);

	// double true_cross_row_similarity = csr_cross_row_similarity(csr->row_ptr, csr->col_ind, csr->nr_rows, csr->nr_cols, csr->nr_nzeros, window_size);

	long num_clusters = csr_clusters_number(csr->row_ptr, csr->col_ind, csr->nr_rows, csr->nr_cols, csr->nr_nzeros, 0);

	// plot_csr(csr, matrix_name);

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
	printf("num_clusters=%ld, ", num_clusters);
	printf("avg_cluster_size=%g, ", ((double) csr->nr_nzeros) / num_clusters);
	printf("\n");

	// fprintf(stderr, "%s\t%ld\t%ld\t%lf\t%lf\t%s\t%s\t%lf\t%lf\n", matrix_name, nr_rows, nr_cols, avg_nnz_per_row, std_nnz_per_row, distribution, placement, bw, skew);
	// fprintf(stderr, "%s\t%d\t%d\t%lf\t%lf\t%s\t%s\t%lf\t%lf\t%lf\t%lf\n", matrix_name, csr->nr_rows, csr->nr_cols, csr->avg_nnz_per_row, csr->std_nnz_per_row, csr->distribution, csr->placement, csr->avg_bw_scaled, csr->skew, csr->avg_num_neighbours, csr->cross_row_similarity);

	fprintf(stderr, "%ld\t%ld\t%lf\t%lf\t%s\t%s\t%lf\t%lf\t%lf\t%lf", nr_rows, nr_cols, avg_nnz_per_row, std_nnz_per_row, distribution, placement, bw, skew, avg_num_neighbours, cross_row_similarity);
	fprintf(stderr, "\t%d\t%d\t%lf\t%lf\t%s\t%s\t%lf\t%lf\t%lf\t%lf\n", csr->nr_rows, csr->nr_cols, csr->avg_nnz_per_row, csr->std_nnz_per_row, csr->distribution, csr->placement, csr->avg_bw_scaled, csr->skew, csr->avg_num_neighbours, csr->cross_row_similarity);

	// fprintf(stderr, "%lf\t%lf\n", csr->avg_num_neighbours, std_neigh);
	// fprintf(stderr, "%lf\n", csr->cross_row_similarity);
	// fprintf(stderr, "%s\t%lf\t%lf\n", matrix_name, csr->avg_num_neighbours/avg_num_neighbours, csr->cross_row_similarity / cross_row_similarity);
	// fprintf(stderr, "%lf\t%lf\n", csr->avg_num_neighbours, csr->cross_row_similarity);
	// fprintf(stderr, "%s\t%lf\t%ld\t%lf\t%lf\n", matrix_name, avg_num_neighbours, num_clusters, 1 / (1 - avg_num_neighbours/2), ((double) csr->nr_nzeros) / num_clusters);
 
	// csr_matrix_write_mtx(csr, "out.mtx");

	free_csr_matrix(csr);
	return 0;
}


