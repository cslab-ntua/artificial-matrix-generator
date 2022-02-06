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


//==========================================================================================================================================
//= Calculate Bandwidth
//==========================================================================================================================================


/*
 * <-----------B----------->
 *         <---x--->
 * [ |...| |O|...|O| |...| ]
 * Possibility of 'n' elements forming an interval of length 'x' in an interval of length 'B' (B > x).
 * a. select the x interval in B                    : B-x+1
 * b. select the 2 elements at the bounds of x      : n*(n-1)
 * c. select the positions of the rest elements     : (x-2)*(x-3)*...*(x-2-(n-2+1))
 * d. all possible positions of the n elements in B : B*(B-1)*...*(B-n+1)
 *
 * pdf(x, B, n) = a*b*c/d = (B-x+1) * (n*(n-1))/(x*(x-1)) * Prod{(x-i)/(B-i), i in [0,n-1]}
 *
 * recursive calculation of product:
 *     f(x) = Prod{(x-i)/(B-i), i in [0,n-1]}
 *     f(x+1) = ((x+1)*(x-n+1)) * f(x)
 */


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
	n = floor(n + 0.5);
	if (n < 2)     // Bandwidth can't be defined with n < 2.
		n = 2;
	E = expected_bw(B, n);
	if (E == 0)     // Bandwidth given was smaller than n.
		B_new = n;
	else
	{
		ratio = B / E;
		B_new = B * ratio;
	}
	debug("n = %g, bw = %g, predicted bw = %g, corrected bw = %g, new prediction = %g\n", n, B, E, B_new, expected_bw(B_new, n));
	return B_new;
}


//==========================================================================================================================================
//= Generator
//==========================================================================================================================================


/*
 * 'bw_scaled'              : Scaled target bandwidth ([0,1]), as a fraction of row size (number of columns).
 * 'avg_num_neighbours'     : Average number of nnz neighbours, i.e. for each nnz, the number of nnz in the left and right adjacent columns ([0, 2]).
 * 'cross_row_similarity'   : Probability of having similar neighbouring rows ([0,1]).
 */
struct csr_matrix *
artificial_matrix_generation(long nr_rows, long nr_cols, double avg_nnz_per_row, double std_nnz_per_row, char * distribution, unsigned int seed, char * placement, double bw_scaled, double skew, double avg_num_neighbours, double cross_row_similarity)
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
	double avg_distr, std_distr;

	double p_right_neigh = avg_num_neighbours / 2;              // Probability of having a right neighbour.
	double expected_cluster_size = 1 / (1 - p_right_neigh);     // Average cluster (adjacent nnz) size, calculated from p_right_neigh (expected value of geometric distribution).

	long t_max_degree[num_threads];
	long max_degree = 0;

	double * degrees;
	double * bandwidths;
	double * scatters;

	debug("arguments: %ld %ld %g %g %s %u %s %g %g %g %g\n", nr_rows, nr_cols, avg_nnz_per_row, std_nnz_per_row, distribution, seed, placement, bw_scaled, skew, avg_num_neighbours, cross_row_similarity);

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
	c_bound1 = MAX / avg_nnz_per_row;                                                          // avg_distr > 0
	c_bound2 = (std_nnz_per_row > 0) ? MAX*MAX / (2 * std_nnz_per_row*std_nnz_per_row) : 0;    // std_distr > 0
	C = (c_bound1 > c_bound2) ? c_bound1 : c_bound2;

	if (C <= 2 || skew < 1)  // For C<=2 the models don't work. Also for skew < 1 we have terrible estimations.
	{
		MAX = 0;
		C = 0;
		avg_exp = 0;
		std_exp = 0;
		avg_distr = avg_nnz_per_row;
		std_distr = std_nnz_per_row;
	}
	else
	{
		avg_exp = MAX / C;
		std_exp = sqrt(avg_exp*avg_exp * (C/2 - 1));

		avg_distr = avg_nnz_per_row - MAX / C;
		std_distr = sqrt(std_nnz_per_row*std_nnz_per_row - std_exp*std_exp);

		if (std_distr > avg_distr / 3)       // For normal distribution, in order to have very few negative random numbers.
			std_distr = avg_distr / 3;
	}

	debug("C = %g\n", C);
	debug("avg exp = %g\n", avg_exp);
	debug("std exp = %g\n", std_exp);
	debug("avg normal = %g\n", avg_distr);
	debug("std normal = %g\n", std_distr);
	debug("halflife = %g\n", log(2) * (nr_rows / C));


	_Pragma("omp parallel")
	{
		int tnum = omp_get_thread_num();
		struct Random_State * rs;
		long i, j, k, i_s, i_e, j_s, j_e, per_t_len;
		double degree_exp, degree_distr;
		double degree;
		long sum, total_sum;
		long local_max_degree = 0;
		struct ordered_set * OS;

		// We anchor the thread iterations, so that the results are independent of the number of threads.
		// 'anchors_num' should be enough for a balanced distribution of rows to threads.
		long anchors_num;
		long anchors_period;
		anchors_num = nr_rows / 100;  // We would like at least 100 rows per part.
		if (anchors_num < 0)
			anchors_num = 1;
		else if (anchors_num > 1<<12) // Enough parts to uniformly distribute to the threads.
			anchors_num = 1<<12;
		anchors_period = nr_rows / anchors_num;

		per_t_len = nr_rows / num_threads;
		i_s = per_t_len * tnum;
		i_s = i_s - i_s % anchors_period;
		if (tnum == num_threads - 1)
			i_e = nr_rows;
		else
		{
			i_e = per_t_len * (tnum + 1);
			i_e = i_e - i_e % anchors_period;
		}

		rs = random_new(tnum);

		// Calculate number of non-zeros for each row.
		sum = 0;
		for (i=i_s;i<i_e;i++)
		{
			if (i % anchors_period == 0)            // Periodic reseeding.
				random_reseed(rs, seed + i);    // We need reproducible results, independently of the number of threads, but random_r() is too slow (even x10 for many, small rows)!
			degree_exp = MAX * exp(-C/nr_rows * (double)i);

			// degree_distr = random_normal(rs, avg_nnz_per_row, std_nnz_per_row);

			if (std_distr > 0)
			{
				degree_distr = random_normal(rs, avg_distr, std_distr);
				// double k, theta;
				// k = (avg_distr*avg_distr) / (std_distr*std_distr);
				// theta = (avg_distr > 0) ? (std_distr*std_distr) / avg_distr : 1;
				// degree_distr = random_gamma(rs, k, theta);
			}
			else
				degree_distr = avg_distr;

			degree = degree_exp + degree_distr;
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
			// bw_corrected = calculate_new_bw(bw_scaled * nr_cols, csr->avg_nnz_per_row);      // Recalculate banwidth with the actual avg nnz per row.
			bw_corrected = calculate_new_bw(bw_scaled * nr_cols, csr->avg_nnz_per_row / expected_cluster_size);      // Recalculate banwidth with the actual avg_nnz_per_row and the expected_cluster_size.

			for (j=1;j<csr->mem_footprint;)
				j <<= 1;
			snprintf(csr->mem_range, sizeof(csr->mem_range), "[%ld-%ld]", j>>1, j);
		}
		_Pragma("omp barrier")

		long bound_l, bound_r, bound_relaxed_l, bound_relaxed_r, range, half_range;         // Relaxed bounds can extend outside the matrix bounds, to spread the average bandwidth if needed.
		double b, s;
		long retries;
		long d1, d2;

		bound_l = 0;
		bound_r = nr_cols;
		OS = ordered_set_new(max_degree);
		j_s = t_base[tnum];
		for (i=i_s;i<i_e;i++)
		{
			if (i % anchors_period == 0)             // Periodic reseeding.
				random_reseed(rs, seed + i);

			degree = degrees[i];
			j_e = j_s + degree;
			offsets[i] = j_s;
			bandwidths[i] = 0;
			scatters[i] = 0;
			if (degree == 0)
				continue;
			OS->size = 0;

			range = floor(bw_corrected + 0.5);
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

			j = j_s;

			if ((cross_row_similarity > 0) && (i % anchors_period != 0) && (random_uniform(rs, 0, 1) <= cross_row_similarity))       // Row similarity.
			{
				long j_prev_s, j_prev_e, degree_prev, degree_min;
				degree_prev = degrees[i-1];
				degree_min = degree < degree_prev ? degree : degree_prev;
				j_prev_s = j_s - degree_prev;
				j_prev_e = j_prev_s + degree_min;
				while (j_prev_s < j_prev_e)
				{
					k = col_ind[j_prev_s];
					// k = col_ind[j_prev_s] + 1;
					if (k >= bound_r)
						k = bound_l;
					if (!ordered_set_insert(OS, k))
					{
						printf("error");
						exit(1);
					}
					values[j] = random_uniform(rs, 0, 1);
					j++;
					j_prev_s++;
				}
			}

			d1 = 0;
			d2 = 0;
			for (;j<j_e;)
			{
				retries = 0;
				do {
					if (retries < 20)   // If we fail 20 times, we can assume it is nearly filled (e.g. 1/2 ^ 20 = 1/1048576 for half filled bandwidth).
					{
						// Random retries give surprisingly better results, compared to searching serially (e.g. 5 vs 45 sec for filled bandwidth), guess better statistical properties.
						k = random_uniform_integer(rs, bound_relaxed_l, bound_relaxed_r);
						if (k >= bound_r)
							k = bound_r - 1;
						else if (k < bound_l)
							k = bound_l;
						retries++;
					}
					else    // Try the bounds, storing the positions to keep O(n) total.
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
				} while (!ordered_set_insert(OS, k));

				values[j] = random_uniform(rs, 0, 1);
				j++;

				if (p_right_neigh > 0)        // Clustering of non-zeros.
				{
					while (j < j_e)
					{
						if (random_uniform(rs, 0, 1) > p_right_neigh)
							break;
						k++;
						if (k >= bound_r)
							k = bound_l;
						if (ordered_set_insert(OS, k))        // Don't retry if the position is taken, because searching for next empty spot can become O(n) hard.
						{
							values[j] = random_uniform(rs, 0, 1);
							j++;
						}
					}
				}
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

	csr->avg_num_neighbours = csr_avg_row_neighbours(csr->row_ptr, csr->col_ind, csr->nr_rows, csr->nr_cols, csr->nr_nzeros, 1);
	csr->cross_row_similarity = csr_cross_row_similarity(csr->row_ptr, csr->col_ind, csr->nr_rows, csr->nr_cols, csr->nr_nzeros, 1);

	free(degrees);
	free(bandwidths);
	free(scatters);

	return csr;
}


//==========================================================================================================================================
//= Utilities
//==========================================================================================================================================


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

