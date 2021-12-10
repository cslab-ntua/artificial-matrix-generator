#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

#include "random.h"
#include "matrix_util.h"
#include "ordered_set.h"

#include "artificial_matrix_generation.h"


#define error(fmt, ...)                       \
do {                                          \
	fprintf(stderr, fmt, __VA_ARGS__);    \
	exit(1);                              \
} while (0)


#ifdef VERBOSE

	#define debug(fmt, ...)              \
	do {                                 \
		printf(fmt, __VA_ARGS__);    \
	} while (0)

#else
	#define debug(...)
#endif


static double mb_list[] =  {4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096};
static long mb_list_n = sizeof(mb_list) / sizeof(mb_list[0]); 


int
free_csr_matrix(struct csr_matrix * csr)
{
	if (!csr)
		return 1;
	free(csr->values);
	free(csr->col_ind);
	free(csr->row_ptr);
	free(csr);
	return 0;
}


static
double
pdf(double x, double B, double n)
{
	double a, b;
	long i;
	a = (B - x + 1) * n/x * (n-1)/(x-1);
	b = 1;
	for (i=0;i<(long)n;i++)
		b *= (x-i) / (B-i);
	return a * b;
}

__attribute__((unused))
static
double
expected_bw_slow(double B, double n)
{
	double E;
	long i;
	E = 0;
	for (i=n;i<=(long)B;i++)
		E += pdf(i, B, n) * (double)i;
	return E;
}


static
double
expected_bw(double B, double n)
{
	double E;
	long i;
	double a, b, b_e, x;

	b_e = 0;
	for (i=0;i<(long)n;i++)
		b_e += log2(n-i) - log2(B-i);
	b = pow(2, b_e);

	E = 0;
	for (i=n;i<=(long)B;i++)
	{
		x = i;
		a = (B - x + 1) * n/x * (n-1)/(x-1);
		E += a * b * (double)i;
		b_e += log2(x+1) - log2(x+1-n);
		b = pow(2, b_e);
	}
	return E;
}


static
double
calculate_new_bw(double B, double n)
{
	double E, ratio, B_new;
	// E = expected_bw_slow(B, n);
	E = expected_bw(B, n);
	if (E == 0)     // Bandwidth given was smaller than n.
		B_new = n;
	else
	{
		ratio = E / B;
		B_new = B / ratio;
	}
	debug("n = %g, bw = %g, predicted bw = %g, corrected bw = %g, new prediction = %g\n", n, B, E, B_new, expected_bw(B_new, n));
	return B_new;
}


/*
 * bw_scaled: scaled target bandwidth ([0,1]), as a fraction of row size (number of columns).
 */
struct csr_matrix *
artificial_matrix_generation(long nr_rows, long nr_cols, double avg_nnz_per_row, double std_nnz_per_row, char * distribution, unsigned int seed, char * placement, double bw_scaled, double skew)
{
	int num_threads = omp_get_max_threads();
	int * offsets;
	int * col_ind;
	long nnz;
	long t_base[num_threads];
	struct csr_matrix * csr;
	ValueType * values;
	double bw_corrected;

	double MAX;
	double C, c_bound1, c_bound2;
	double avg_exp, std_exp;
	double avg_norm, std_norm;

	long t_max_degree[num_threads];
	long max_degree = 0;

	double * degrees;
	double * bandwidths;
	double * scatters;

	debug("arguments: %ld %ld %g %g %s %u %s %g %g\n", nr_rows, nr_cols, avg_nnz_per_row, std_nnz_per_row, distribution, seed, placement, bw_scaled, skew);

	offsets = (typeof(offsets)) malloc((nr_rows+1) * sizeof(*offsets));
	csr = (typeof(csr)) malloc(sizeof(*csr));
	csr->nr_rows = nr_rows;
	csr->nr_cols = nr_cols;
	csr->seed = seed;
	csr->distribution = distribution;
	csr->placement = placement;
	if (bw_scaled > 1)                    // bandwidth <= number of columns
		bw_scaled = 1;

	degrees = (typeof(degrees)) malloc(nr_rows * sizeof(*degrees));
	bandwidths = (typeof(bandwidths)) malloc(nr_rows * sizeof(*bandwidths));
	scatters = (typeof(scatters)) malloc(nr_rows * sizeof(*scatters));


	// Calculate parameters from the addition of the exponential.
	MAX = avg_nnz_per_row * (1 + skew);
	c_bound1 = MAX / avg_nnz_per_row;
	c_bound2 = (std_nnz_per_row > 0) ? MAX*MAX / (2 * std_nnz_per_row*std_nnz_per_row) : 0;
	C = (c_bound1 > c_bound2) ? c_bound1 : c_bound2;

	if (C <= 2 || skew < 1)  // For C<=2 the models don't work. Also for skew < 1 we have terrible estimations.
	{
		MAX = 0;
		C = 0;
		avg_exp = 0;
		std_exp = 0;
		avg_norm = avg_nnz_per_row;
		std_norm = std_nnz_per_row;
	}
	else
	{
		avg_exp = MAX / C;
		std_exp = sqrt(avg_exp*avg_exp * (C/2 - 1));

		avg_norm = avg_nnz_per_row - MAX / C;
		std_norm = sqrt(std_nnz_per_row*std_nnz_per_row - std_exp*std_exp);
		if (avg_norm / 3 < std_norm)
			std_norm = avg_norm / 3;
	}

	debug("C = %g\n", C);
	debug("avg exp = %g\n", avg_exp);
	debug("std exp = %g\n", std_exp);
	debug("avg normal = %g\n", avg_norm);
	debug("std normal = %g\n", std_norm);
	debug("halflife = %g\n", log(2) * (nr_rows / C));


	_Pragma("omp parallel")
	{
		int tnum = omp_get_thread_num();
		struct Random_State * rs;
		long reseed_period;
		long i, j, k, i_s, i_e, j_s, j_e, per_t_len;
		double e, norm;
		double degree;
		long sum, total_sum;
		long local_max_degree = 0;
		struct ordered_set * OS;

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
			e = MAX * exp(-C/nr_rows * (double)i);
			// norm = random_normal(rs, avg_nnz_per_row, std_nnz_per_row);
			norm = random_normal(rs, avg_norm, std_norm);
			degree = e + norm;
			degree = (degree > 0) ? floor(degree + 0.5) : 0;
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
			long tmp;
			total_sum = 0;
			for (i=0;i<num_threads;i++)
			{
				tmp = t_base[i];
				t_base[i] = total_sum;
				total_sum += tmp;
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

			csr->avg_nnz_per_row = ((double) nnz) / nr_rows;
			bw_corrected = calculate_new_bw(bw_scaled * nr_cols, floor(csr->avg_nnz_per_row + 0.5));      // Recalculate banwidth with the actual avg nnz per row.

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

		long bound_l, bound_r, bound_relaxed_l, bound_relaxed_r, range, half_range;
		double b, s;
		long retries;
		long d1, d2;

		bound_l = 0;
		bound_r = nr_cols;
		OS = ordered_set_new(max_degree);
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
			OS->size = 0;

			range = floor(bw_corrected + 0.5);
			// range = floor(bw_scaled * nr_cols + 0.5);
			if (range < degree)                           // At this point: range >= degree > 0
				range = degree;
			half_range = range / 2;

			if (*placement == 'd')
			{
				bound_relaxed_l = i - half_range;
				bound_relaxed_r = i + half_range + range % 2;      // Correct modulo upwards.
			}
			else
			{
				bound_relaxed_l = nr_cols/2 - half_range;
				bound_relaxed_r = nr_cols/2 + half_range + range % 2;
			}
			bound_l = bound_relaxed_l;
			bound_r = bound_relaxed_r;
			if (bound_l < 0)
			{
				bound_l = 0;
				if (bound_r < degree)
					bound_r = degree;
			}
			if (bound_r > nr_cols)
			{
				bound_r = nr_cols;
				if (bound_l > nr_cols - degree)
					bound_l = nr_cols - degree;
			}

			if (i % reseed_period == 0)
				random_reseed(rs, seed + i);
			d1 = 0;
			d2 = 0;
			for (j=j_s;j<j_e;j++)
			{
				k = random_uniform_integer(rs, bound_relaxed_l, bound_relaxed_r);
				if (k >= bound_r)
					k = bound_r - 1;
				else if (k < bound_l)
					k = bound_l;
				retries = 0;
				while (!ordered_set_insert(OS, k))
				{
					retries++;
					if (retries < 20)   // If we fail 20 times, we can assume it is nearly filled (e.g. 1/2 ^ 20 = 1/1048576 for half filled bandwidth).
					{
						// Random retries give surprisingly better results, compared to searching serially (e.g. 5 vs 45 sec for filled bandwidth), guess better statistical properties.
						k = random_uniform_integer(rs, bound_relaxed_l, bound_relaxed_r);
						if (k >= bound_r)
							k = bound_r - 1;
						else if (k < bound_l)
							k = bound_l;
						// k = random_uniform_integer(rs, bound_l, bound_r);
						// k++;
						// if (k >= bound_r)
							// k = bound_l;
					}
					else
					{
						if (d1 < d2)
						{
							k = bound_l + d1;
							d1++;
						}
						else
						{
							k = bound_r - d2 - 1;
							d2++;
						}
					}
				}

				values[j] = random_uniform(rs, 0, 1);
			}

			ordered_set_sort(OS, &col_ind[j_s]);

			b = col_ind[j_e-1] - col_ind[j_s] + 1;
			bandwidths[i] = b;

			// s = (b > 0) ? degree / b : 0;
			s = degree / b;
			scatters[i] = s;

			j_s = j_e;
		}

		ordered_set_destroy(OS);
		random_destroy(rs);
	}

	csr->avg_nnz_per_row = ((double) nnz) / nr_rows;
	csr->std_nnz_per_row = matrix_std_base(degrees, nr_rows, csr->avg_nnz_per_row);
	matrix_min_max(degrees, nr_rows, &csr->min_nnz_per_row, &csr->max_nnz_per_row);
	csr->skew = (csr->max_nnz_per_row - csr->avg_nnz_per_row) / csr->avg_nnz_per_row;

	csr->avg_bw = matrix_mean(bandwidths, nr_rows);
	csr->std_bw = matrix_std_base(bandwidths, nr_rows, csr->avg_bw);
	csr->avg_bw_scaled = csr->avg_bw / nr_cols;
	csr->std_bw_scaled = csr->std_bw / nr_cols;

	csr->avg_sc = matrix_mean(scatters, nr_rows);
	csr->std_sc = matrix_std_base(scatters, nr_rows, csr->avg_sc);
	csr->avg_sc_scaled = csr->avg_sc * nr_cols;
	csr->std_sc_scaled = csr->std_sc * nr_cols;

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

