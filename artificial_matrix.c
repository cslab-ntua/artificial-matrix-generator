#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

#include "artificial_matrix_generation.h"
#include "matrix_util.h"


#include "time_it.h"
#include "plot/plot.h"

double
get_row(void * A, int pos)
{
	struct csr_matrix * csr = A;
	long i = binary_search(csr->row_ptr, 0, csr->nr_rows, pos, NULL, NULL);
	if (csr->row_ptr[i] > pos)
		i--;
	return (double) i;
}


double
get_col(void * A, int pos)
{
	struct csr_matrix * csr = A;
	return (double) csr->col_ind[pos];
}


double
get_degree(void * A, int i)
{
	struct csr_matrix * csr = A;
	return (double) csr->row_ptr[i+1] -  csr->row_ptr[i];
}


void
plot_csr(struct csr_matrix * csr, char * matrix_name, double * neigh_distances_freq_perc)
{
	long num_pixels = 1024;
	long num_pixels_x = num_pixels, num_pixels_y = num_pixels;
	long buf_n = strlen(matrix_name) + 1 + 1000;
	char buf[buf_n], buf_title[buf_n];

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
	long num_bins = 10000;
	snprintf(buf, buf_n, "figures/%s-artificial_degree_distribution.png", matrix_name);
	snprintf(buf_title, buf_n, "%s-artificial: degree distribution (percentages)", matrix_name);
	figure_simple_plot(buf, num_pixels_x, num_pixels_y, (NULL, csr, NULL, csr->nr_rows, 0, NULL, get_degree),
		figure_enable_legend(_fig);
		figure_set_title(_fig, buf_title);
		figure_series_type_histogram(_s, num_bins, 1);
		figure_series_type_barplot(_s);
	);

	// Plot neighbour distances frequencies.
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
}


double *
csr_row_neighbours(struct csr_matrix * csr, long window_size)
{
	double * neigh_num = malloc(csr->nr_nzeros * sizeof(*neigh_num));
	#pragma omp parallel
	{
		long i, j, k;
		#pragma omp for schedule(static)
		for (i=0;i<csr->nr_nzeros;i++)
			neigh_num[i] = 0;
		#pragma omp for schedule(static)
		for (i=0;i<csr->nr_rows;i++)
		{
			for (j=csr->row_ptr[i];j<csr->row_ptr[i+1];j++)
			{
				for (k=j+1;k<csr->row_ptr[i+1];k++)
				{
					if (csr->col_ind[k] - csr->col_ind[j] > window_size)
						break;
					neigh_num[j]++;
					neigh_num[k]++;
				}
			}
		}
	}
	return neigh_num;
}


double *
csr_neighbours_distances_frequencies(struct csr_matrix * csr)
{
	long * frequencies = malloc(csr->nr_cols * sizeof(*frequencies));
	double * frequencies_percentage = malloc(csr->nr_cols * sizeof(*frequencies_percentage));
	#pragma omp parallel
	{
		long i, j;
		#pragma omp for schedule(static)
		for (i=0;i<csr->nr_cols;i++)
			frequencies[i] = 0;
		#pragma omp for schedule(static)
		for (i=0;i<csr->nr_rows;i++)
		{
			for (j=csr->row_ptr[i];j<csr->row_ptr[i+1]-1;j++)
			{
				__atomic_fetch_add(&frequencies[csr->col_ind[j+1] - csr->col_ind[j]], 1, __ATOMIC_RELAXED);
			}
		}
		#pragma omp for schedule(static)
		for (i=0;i<csr->nr_cols;i++)
			frequencies_percentage[i] = ((double) 100 * frequencies[i]) / csr->nr_nzeros;
	}
	free(frequencies);
	return frequencies_percentage;
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
		csr = artificial_matrix_generation(nr_rows, nr_cols, avg_nnz_per_row, std_nnz_per_row, distribution, seed, placement, bw, skew, avg_num_neighbours);
	);
	printf("time generate matrix = %g\n", time);

	long window_size;  // Distance from left and right.
	// window_size = 64 / sizeof(double) / 2;
	// window_size = 8;
	// window_size = 4;
	window_size = 1;
	double * neigh_num = csr_row_neighbours(csr, window_size);
	double mean_neigh = matrix_mean(neigh_num, csr->nr_nzeros);
	double std_neigh = matrix_std_base(neigh_num, csr->nr_nzeros, mean_neigh);

	// double * neigh_distances_freq_perc = csr_neighbours_distances_frequencies(csr);
	// plot_csr(csr, matrix_name, neigh_distances_freq_perc);

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
	printf("mean neigh num = %g, ", mean_neigh);
	printf("std neigh num = %g, ", std_neigh);
	printf("\n");

	// fprintf(stderr, "%s\t%ld\t%ld\t%lf\t%lf\t%s\t%s\t%lf\t%lf\n", matrix_name, nr_rows, nr_cols, avg_nnz_per_row, std_nnz_per_row, distribution, placement, bw, skew);
	// fprintf(stderr, "\t%d\t%d\t%lf\t%lf\t%s\t%s\t%lf\t%lf\n", csr->nr_rows, csr->nr_cols, csr->avg_nnz_per_row, csr->std_nnz_per_row, csr->distribution, csr->placement, csr->avg_bw_scaled, csr->skew);

	// fprintf(stderr, "%s\t%lf\t%lf\n", matrix_name, mean_neigh, std_neigh);
	fprintf(stderr, "%lf\t%lf\n", mean_neigh, std_neigh);
 
	// csr_matrix_write_mtx(csr, "out.mtx");

	free_csr_matrix(csr);
	return 0;
}


