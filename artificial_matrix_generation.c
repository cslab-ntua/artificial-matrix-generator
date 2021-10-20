#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

#include "random.h"
#include "matrix_util.h"
#include "sorted_set.h"

#include "artificial_matrix_generation.h"


struct csr_matrix *
artificial_matrix_generation(long nr_rows, long nr_cols, double avg_nnz_per_row, double std_nnz_per_row, char * distribution, unsigned int seed, char * placement, double d_f)
{
	int num_threads = omp_get_max_threads();
	int * offsets;
	int * col_ind;
	long nnz;
	long t_base[num_threads];
	struct csr_matrix * csr;

	long t_max_degree[num_threads];
	long max_degree = 0;

	double * degrees;
	double * bandwidths;
	double * scatters;

	long i;

	offsets = malloc((nr_rows+1) * sizeof(*offsets));
	csr = malloc(sizeof(*csr));
	csr->nr_rows = nr_rows;
	csr->nr_cols = nr_cols;
	csr->seed = seed;
	csr->placement = placement;
	csr->diagonal_factor = d_f;

	degrees = malloc(nr_rows * sizeof(*degrees));
	bandwidths = malloc(nr_rows * sizeof(*bandwidths));
	scatters = malloc(nr_rows * sizeof(*scatters));

	_Pragma("omp parallel")
	{
		int tnum = omp_get_thread_num();
		struct Random_State * rs;
		long reseed_period;
		long i, j, k, i_s, i_e, j_s, j_e, per_t_len;
		double rand_val;
		long degree;
		long sum, total_sum;
		long local_max_degree = 0;
		struct sorted_set * SU;

		reseed_period = nr_rows / 1000;
		if (reseed_period < 1)
			reseed_period = 1;

		per_t_len = nr_rows / num_threads;
		i_s = per_t_len * tnum;
		i_s = i_s - i_s % reseed_period;
		if (tnum == num_threads - 1)
			i_e = nr_rows;
		else
		{
			i_e = per_t_len * (tnum + 1);
			i_e = i_e - i_e % reseed_period;
		}

		rs = random_new(tnum);

		sum = 0;
		for (i=i_s;i<i_e;i++)
		{
			if (i % reseed_period == 0)
				random_reseed(rs, seed + i);    // We need reproducible results, independently of the number of threads, but random_r() is too slow (even x10 for many, small rows)!
			rand_val = random_normal(rs, avg_nnz_per_row, std_nnz_per_row);
			degree = (long) (rand_val > 0 ? rand_val : -rand_val);
			if (degree > nr_cols)
				degree = nr_cols;
			if (degree > local_max_degree)
				local_max_degree = degree;
			degrees[i] = degree;
			offsets[i] = degree;
			sum += degree;
		}

		__atomic_store_n(&t_base[tnum], sum, __ATOMIC_RELAXED);
		__atomic_store_n(&t_max_degree[tnum], local_max_degree, __ATOMIC_RELAXED);

		_Pragma("omp barrier")
		_Pragma("omp single")
		{
			total_sum = 0;
			for (i=0;i<num_threads;i++)
			{
				degree = t_base[i];
				t_base[i] = total_sum;
				total_sum += degree;

				if (t_max_degree[i] > max_degree)
					max_degree = t_max_degree[i];
			}
			nnz = total_sum;
			offsets[nr_rows] = nnz;
			col_ind = malloc(nnz * sizeof(*col_ind));

			csr->row_ptr = offsets;
			csr->col_ind = col_ind;
			csr->nr_nzeros = nnz;
		}
		_Pragma("omp barrier")

		long bound_l, bound_r, half_range;
		double b, s;
		bound_l = 0;
		bound_r = nr_cols;
		SU = sorted_set_new(max_degree);
		j_s = t_base[tnum];
		for (i=i_s;i<i_e;i++)
		{
			degree = offsets[i];
			j_e = j_s + degree;
			offsets[i] = j_s;
			bandwidths[i] = 0;
			scatters[i] = 0;
			if (degree == 0)
				continue;
			SU->size = 0;
			if (*placement == 'd')
			{
				half_range = (((double) degree) / d_f);
				bound_l = i - half_range;
				if (bound_l < 0)
					bound_l = 0;
				bound_r = i + half_range;
				if (bound_r > nr_cols)
					bound_r = nr_cols;
			}
			if (i % reseed_period == 0)
				random_reseed(rs, seed + i);
			// long tmp = 0;
			for (j=j_s;j<j_e;j++)
			{
				k = random_uniform_integer(rs, bound_l, bound_r);
				while (!sorted_set_insert(SU, k))
				{
					// k = random_uniform_integer(rs, bound_l, bound_r);
					k++;
					if (k >= bound_r)
						k = bound_l;
					// tmp++;
				}
			}
			// if (tmp > 10)
				// printf("tmp=%ld\n", tmp);

			sorted_set_sort(SU, &col_ind[j_s]);

			b = col_ind[j_e-1] - col_ind[j_s];
			bandwidths[i] = b;

			s = (b > 0) ? degree / b : 0;
			scatters[i] = s;

			j_s = j_e;
		}

		sorted_set_destroy(SU);
		random_destroy(rs);
	}

	csr->avg_nnz_per_row = ((double) nnz) / nr_rows;
	csr->std_nnz_per_row = matrix_std_base(degrees, nr_rows, csr->avg_nnz_per_row);

	csr->avg_bw = matrix_mean(bandwidths, nr_rows);
	csr->std_bw = matrix_std_base(bandwidths, nr_rows, csr->avg_bw);

	csr->avg_sc = matrix_mean(scatters, nr_rows);
	csr->std_sc = matrix_std_base(scatters, nr_rows, csr->avg_sc);

	for (i=0;i<10;i++)
		printf("%d ", offsets[i]);
	printf("\n");

	#if 0
		printf("%d %d %d %d\n", csr->nr_rows, csr->nr_cols, csr->nr_nzeros, csr->row_ptr[nr_rows]);
		for (i=0;i<csr->nr_rows;i++)
		{
			long j, j_s, j_e;
			j_s = csr->row_ptr[i];
			j_e = csr->row_ptr[i+1];
			printf("%3ld (%3ld): ", i, j_e - j_s);
			for (j=j_s;j<j_e;j++)
			{
				printf("%3d ", csr->col_ind[j]);
			}
			printf("\n");
		}
	#endif

	free(degrees);
	free(bandwidths);
	free(scatters);

	return csr;
}

