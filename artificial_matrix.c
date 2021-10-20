#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

#include "time_it.h"

#include "artificial_matrix_generation.h"


#ifdef PLOT

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
		);

		// Plot degree histogram.
		#if 1
			long resolution = 1000;
			figure_simple_plot("degree_distribution.png", num_pixels_x, num_pixels_y, (NULL, csr, NULL, csr->nr_rows, 0, NULL, get_degree),
				figure_enable_legend(_fig);
				figure_set_title(_fig, "degree distribution");
				figure_series_type_histogram(_s, resolution);
				// figure_set_bounds_x(_fig, -3, 5);
			);
		#endif
	}

#endif


int
main()
{
	struct csr_matrix * csr;
	long nr_rows = 1000000;
	long nr_cols = nr_rows;
	double avg_nnz_per_row, std_nnz_per_row;
	unsigned int seed = 0;
	char * placement;
	double d_f = 0.05;
	double time;

	// placement = "random";
	placement = "diagonal";

	avg_nnz_per_row = 10;
	std_nnz_per_row = 2.5;
	time = time_it(1,
		csr = artificial_matrix_generation(nr_rows, nr_cols, avg_nnz_per_row, std_nnz_per_row, seed, placement, d_f);
	);
	printf("time generate matrix = %g\n", time);

	#ifdef PLOT
		plot_csr(csr, "csr.png");
	#endif

	printf("rows=%d, ", csr->nr_rows);
	printf("cols=%d, ", csr->nr_cols);
	printf("nnz=%d, ", csr->nr_nzeros);
	printf("avg_nnz_per_row=%g, ", csr->avg_nnz_per_row);
	printf("std_nnz_per_row=%g, ", csr->std_nnz_per_row);
	printf("avg_bw=%g, ", csr->avg_bw);
	printf("std_bw=%g, ", csr->std_bw);
	printf("avg_sc=%g, ", csr->avg_sc);
	printf("std_sc=%g, ", csr->std_sc);
	printf("\n");
 

	return 0;
}


