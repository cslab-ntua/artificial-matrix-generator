#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

#include "artificial_matrix_generation.h"


#ifdef PLOT

	#include "time_it.h"
	#include "plot/plot.h"

	double
	get_row(void * A, int pos)
	{
		struct csr_matrix * csr = A;
		long i = binary_search(csr->row_ptr, 0, csr->nr_rows, pos);
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
		// long i = binary_search(csr->row_ptr, 0, csr->nr_rows, pos);
		// if (csr->row_ptr[i] > pos)
			// i--;
		return (double) csr->row_ptr[i+1] -  csr->row_ptr[i];
	}


	void
	plot_csr(struct csr_matrix * csr, char * file_out)
	{
		long num_pixels = 1024;
		long num_pixels_x, num_pixels_y;

		if (csr->nr_cols < 1024)
		{
			num_pixels_x = csr->nr_cols;
			num_pixels_y = csr->nr_rows;
		}
		else
		{
			num_pixels_x = num_pixels;
			num_pixels_y = num_pixels;
		}

		figure_simple_plot(file_out, num_pixels_x, num_pixels_y, (csr, csr, NULL, csr->nr_nzeros, 0, get_col, get_row),
			figure_enable_legend(_fig);
			figure_axes_origin_upper_left(_fig);
			figure_set_bounds_x(_fig, 0, csr->nr_cols);
			figure_set_bounds_y(_fig, 0, csr->nr_rows);
		);

		// Plot degree histogram.
		#if 1
			long num_bins = 1000;
			figure_simple_plot("degree_distribution.png", num_pixels_x, num_pixels_y, (NULL, csr, NULL, csr->nr_rows, 0, NULL, get_degree),
				figure_enable_legend(_fig);
				figure_set_title(_fig, "degree distribution");
				figure_series_type_histogram(_s, num_bins, 1);
				// figure_set_bounds_x(_fig, -3, 5);
			);
		#endif
	}

#endif


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
	seed = atoi(argv[i++]);
	if (argc >= i)
	{
		matrix_name = argv[i++];
		printf("matrix: %s\n", matrix_name);
	}


	#ifdef PLOT
	double time;
	time = time_it(1,
	#endif
		csr = artificial_matrix_generation(nr_rows, nr_cols, avg_nnz_per_row, std_nnz_per_row, distribution, seed, placement, bw, skew);
	#ifdef PLOT
	);
	printf("time generate matrix = %g\n", time);
	#endif

	// csr_matrix_write_mtx(csr, "out.mtx");

	#ifdef PLOT
		// plot_csr(csr, "csr.png");
	#endif

	printf("synthetic, ");
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
	printf("\n");

	fprintf(stderr, "%s\t%ld\t%ld\t%lf\t%lf\t%s\t%s\t%lf\t%lf\t%d\n", matrix_name, nr_rows, nr_cols, avg_nnz_per_row, std_nnz_per_row, distribution, placement, bw, skew, seed);
	fprintf(stderr, "\t%d\t%d\t%lf\t%lf\t%s\t%s\t%lf\t%lf\t%d\n", csr->nr_rows, csr->nr_cols, csr->avg_nnz_per_row, csr->std_nnz_per_row, csr->distribution, csr->placement, csr->avg_bw_scaled, csr->skew, csr->seed);
 
	free_csr_matrix(csr);
	return 0;
}


