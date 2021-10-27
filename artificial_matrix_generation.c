#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

#include "random.h"
#include "matrix_util.h"
#include "sorted_set.h"

#include "artificial_matrix_generation.h"


#define error(fmt, ...)                      \
do {                                         \
	printf(stderr, fmt, __VA_ARGS__);    \
	exit(1);                             \
} while (0)


static double mb_list[] =  {4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096};
static long mb_list_n = sizeof(mb_list) / sizeof(mb_list[0]); 


struct csr_matrix *
artificial_matrix_generation(long nr_rows, long nr_cols, double avg_nnz_per_row, double std_nnz_per_row, char * distribution, unsigned int seed, char * placement, double d_f)
{
	int num_threads = omp_get_max_threads();
	int * offsets;
	int * col_ind;
	long nnz;
	long t_base[num_threads];
	struct csr_matrix * csr;
	ValueType * values;

	long t_max_degree[num_threads];
	long max_degree = 0;

	double * degrees;
	double * bandwidths;
	double * scatters;

	offsets = (typeof(offsets)) malloc((nr_rows+1) * sizeof(*offsets));
	csr = (typeof(csr)) malloc(sizeof(*csr));
	csr->nr_rows = nr_rows;
	csr->nr_cols = nr_cols;
	csr->seed = seed;
	csr->distribution = distribution;
	csr->placement = placement;
	csr->diagonal_factor = d_f;

	degrees = (typeof(degrees)) malloc(nr_rows * sizeof(*degrees));
	bandwidths = (typeof(bandwidths)) malloc(nr_rows * sizeof(*bandwidths));
	scatters = (typeof(scatters)) malloc(nr_rows * sizeof(*scatters));

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
		struct sorted_set * SS;

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
			col_ind = (typeof(col_ind)) malloc(nnz * sizeof(*col_ind));
			values = (typeof(values)) malloc(nnz * sizeof(*values));

			csr->nr_nzeros = nnz;
			csr->row_ptr = offsets;
			csr->col_ind = col_ind;
			csr->values = values;
			csr->density = ((double) nnz) / ((double) nr_rows * nr_cols) * 100;
			csr->mem_footprint = (nnz * (sizeof(*values) + sizeof(*col_ind)) +  (nr_rows + 1) * sizeof(*offsets)) / ((double) 1024 * 1024);

			for (i=0;i<mb_list_n;i++)
				if (csr->mem_footprint < mb_list[i])
					break;
			if (i == 0)
				snprintf(csr->mem_range, sizeof(csr->mem_range), "[<%g]", mb_list[0]);
			else if (i >= mb_list_n)
				snprintf(csr->mem_range, sizeof(csr->mem_range), "[>%g]", mb_list[mb_list_n - 1]);
			else
				snprintf(csr->mem_range, sizeof(csr->mem_range), "[%g-%g]", mb_list[i-1], mb_list[i]);
		}
		_Pragma("omp barrier")

		long bound_l, bound_r, half_range;
		double b, s;
		bound_l = 0.0;
		bound_r = nr_cols;
		SS = sorted_set_new(max_degree);
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
			SS->size = 0;
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
				while (!sorted_set_insert(SS, k))
				{
					// k = random_uniform_integer(rs, bound_l, bound_r);
					k++;
					if (k >= bound_r)
						k = bound_l;
					// tmp++;
				}

				values[j] = random_uniform(rs, 0, 1);
			}
			// if (tmp > 10)
				// printf("tmp=%ld\n", tmp);

			sorted_set_sort(SS, &col_ind[j_s]);

			b = col_ind[j_e-1] - col_ind[j_s] + 1;
			b /= nr_cols;
			bandwidths[i] = b;

			// s = (b > 0) ? degree / b : 0;
			s = degree / b;
			scatters[i] = s;

			j_s = j_e;
		}

		sorted_set_destroy(SS);
		random_destroy(rs);
	}

	csr->avg_nnz_per_row = ((double) nnz) / nr_rows;
	csr->std_nnz_per_row = matrix_std_base(degrees, nr_rows, csr->avg_nnz_per_row);

	csr->avg_bw = matrix_mean(bandwidths, nr_rows);
	csr->std_bw = matrix_std_base(bandwidths, nr_rows, csr->avg_bw);

	csr->avg_sc = matrix_mean(scatters, nr_rows);
	csr->std_sc = matrix_std_base(scatters, nr_rows, csr->avg_sc);

	free(degrees);
	free(bandwidths);
	free(scatters);

	return csr;
}


void
csr_matrix_print(struct csr_matrix * csr)
{
	long i, j, j_s, j_e;
	printf("%d %d %d %d\n", csr->nr_rows, csr->nr_cols, csr->nr_nzeros, csr->row_ptr[csr->nr_rows]);
	for (i=0;i<csr->nr_rows;i++)
	{
		j_s = csr->row_ptr[i];
		j_e = csr->row_ptr[i+1];
		printf("%3ld (%3ld): ", i, j_e - j_s);
		for (j=j_s;j<j_e;j++)
		{
			printf("(%3d:%g) ", csr->col_ind[j], csr->values[j]);
		}
		printf("\n");
	}
}


void
csr_matrix_write_mtx(struct csr_matrix * csr, char * file_out)
{
	long i, j, j_s, j_e;
	FILE * file;
	file = fopen(file_out, "w");
	fprintf(file, "%%%%MatrixMarket matrix coordinate real general\n");
	fprintf(file, "%d %d %d\n", csr->nr_rows, csr->nr_cols, csr->nr_nzeros);
	for (i=0;i<csr->nr_rows;i++)
	{
		j_s = csr->row_ptr[i];
		j_e = csr->row_ptr[i+1];
		for (j=j_s;j<j_e;j++)
			fprintf(file, "%ld %ld %g\n", i, j, csr->values[j]);
	}
}

