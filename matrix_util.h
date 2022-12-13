#ifndef MATRIX_UTIL_H
#define MATRIX_UTIL_H

#include <omp.h>
#include <math.h>
#include <stdlib.h>

// #include "artificial_matrix_generation.h"


__attribute__((unused))
static
void
matrix_min_max(double * a, long n, double * min, double * max)
{
	int num_threads = omp_get_max_threads();
	double min_t[num_threads], max_t[num_threads];
	long i;
	#pragma omp parallel
	{
		int tnum = omp_get_thread_num();
		double min_partial, max_partial;
		long i;
		min_partial = a[0];
		max_partial = a[0];
		#pragma omp for schedule(static)
		for (i=0;i<n;i++)
		{
			if (a[i] < min_partial)
				min_partial = a[i];
			if (a[i] > max_partial)
				max_partial = a[i];
		}
		min_t[tnum] = min_partial;
		max_t[tnum] = max_partial;
	}
	*min = min_t[0];
	*max = max_t[0];
	for (i=0;i<num_threads;i++)
	{
		if (min_t[i] < *min)
			*min = min_t[i];
		if (max_t[i] > *max)
			*max = max_t[i];
	}
}


__attribute__((unused))
static
double
matrix_mean(double * a, long n)
{
	int num_threads = omp_get_max_threads();
	double sums[num_threads];
	double sum = 0;
	long i;
	#pragma omp parallel
	{
		int tnum = omp_get_thread_num();
		long i;
		double t_sum = 0;
		#pragma omp for schedule(static)
		for (i=0;i<n;i++)
			t_sum += a[i];
		sums[tnum] = t_sum;
	}
	for (i=0;i<num_threads;i++)
		sum += sums[i];
	return sum / n;
}


__attribute__((unused))
static
double
matrix_var_base(double * a, long n, double mean)
{
	int num_threads = omp_get_max_threads();
	double sums[num_threads];
	double sum = 0;
	long i;
	#pragma omp parallel
	{
		int tnum = omp_get_thread_num();
		long i;
		double t_sum = 0, tmp;
		#pragma omp for schedule(static)
		for (i=0;i<n;i++)
		{
			tmp = a[i] - mean;
			t_sum += tmp * tmp;
		}
		sums[tnum] = t_sum;
	}
	for (i=0;i<num_threads;i++)
		sum += sums[i];
	return sum / n;
}


__attribute__((unused))
static
double
matrix_var(double * a, long n)
{
	double mean = matrix_mean(a, n);
	return matrix_var_base(a, n, mean);
}


__attribute__((unused))
static
double
matrix_std_base(double * a, long n, double mean)
{
	return sqrt(matrix_var_base(a, n, mean));
}


__attribute__((unused))
static
double
matrix_std(double * a, long n)
{
	return matrix_std_base(a, n, matrix_mean(a, n));
}


__attribute__((unused))
static
double *
csr_neighbours_distances_frequencies(int * R_offsets, int * C, long m, long n, long nnz, int ignore_big_rows)
{
	long * frequencies = (typeof(frequencies)) malloc(n * sizeof(*frequencies));
	double * frequencies_percentages = (typeof(frequencies_percentages)) malloc(n * sizeof(*frequencies_percentages));
	double nnz_per_row_avg = (double) nnz / (double) m;
	#pragma omp parallel
	{
		long i, j;
		long degree;
		#pragma omp for schedule(static)
		for (i=0;i<n;i++)
			frequencies[i] = 0;
		#pragma omp for schedule(static)
		for (i=0;i<m;i++)
		{
			degree = R_offsets[i+1] - R_offsets[i];
			if (degree <= 0)
				continue;
			if (ignore_big_rows && (degree > 100 * nnz_per_row_avg))            // Filter out big rows.
			{
				__atomic_fetch_add(&frequencies[0], degree, __ATOMIC_RELAXED);
				continue;
			}
			for (j=R_offsets[i];j<R_offsets[i+1]-1;j++)
			{
				__atomic_fetch_add(&frequencies[C[j+1] - C[j]], 1, __ATOMIC_RELAXED);
			}
			__atomic_fetch_add(&frequencies[0], 1, __ATOMIC_RELAXED);     // Add the last element of the row to the 'no-neighbour' (zero) frequency.
		}
		#pragma omp for schedule(static)
		for (i=0;i<n;i++)
			frequencies_percentages[i] = ((double) 100 * frequencies[i]) / nnz;
	}
	free(frequencies);
	return frequencies_percentages;
}


__attribute__((unused))
static
long
csr_clusters_number(int * R_offsets, int * C, long m, __attribute__((unused)) long n, __attribute__((unused)) long nnz, long max_gap_size)
{
	long total_clusters = 0;
	#pragma omp parallel
	{
		long i, j, k, degree;
		long num_clusters = 0;
		#pragma omp for schedule(static)
		for (i=0;i<m;i++)
		{
			degree = R_offsets[i+1] - R_offsets[i];
			if (degree <= 0)
				continue;
			j = R_offsets[i];
			while (j < R_offsets[i+1])
			{
				k = j + 1;
				while ((k < R_offsets[i+1]) && (C[k] - C[k-1] <= max_gap_size + 1))   // distance 1 means gap 0
					k++;
				num_clusters++;
				j = k;
			}
		}
		__atomic_fetch_add(&total_clusters, num_clusters, __ATOMIC_RELAXED);
	}
	return total_clusters;
}


__attribute__((unused))
static
double
csr_avg_row_neighbours(int * R_offsets, int * C, long m, __attribute__((unused)) long n, long nnz, long window_size)
{
	long total_num_neigh = 0;
	#pragma omp parallel
	{
		long i, j, k;
		long num_neigh = 0;
		#pragma omp for schedule(static)
		for (i=0;i<m;i++)
		{
			for (j=R_offsets[i];j<R_offsets[i+1];j++)
			{
				for (k=j+1;k<R_offsets[i+1];k++)
				{
					if (C[k] - C[j] > window_size)
						break;
					num_neigh += 2;
				}
			}
		}
		__atomic_fetch_add(&total_num_neigh, num_neigh, __ATOMIC_RELAXED);
	}
	return ((double) total_num_neigh) / nnz;
}


__attribute__((unused))
static
void
csr_num_neighbours(int * R_offsets, int * C, long m, __attribute__((unused)) long n, long nnz, long window_size,
	double * avg_num_neighbours, double * std_num_neighbours, double * min_num_neighbours, double * max_num_neighbours)
{
	double * num_neighbours = (typeof(num_neighbours)) malloc(nnz * sizeof(*num_neighbours));
	#pragma omp parallel
	{
		long i, j, k;
		#pragma omp for schedule(static)
		for (i=0;i<m;i++)
		{
			for (j=R_offsets[i];j<R_offsets[i+1];j++)
			{
				long num_neigh_curr = 0;
				for (k=j+1;k<R_offsets[i+1];k++)
				{
					if (C[k] - C[j] > window_size)
						break;
					num_neigh_curr += 2;
				}
				num_neighbours[j] = (double)num_neigh_curr;
			}
		}
	}

	matrix_min_max(num_neighbours, nnz, min_num_neighbours, max_num_neighbours);
	*avg_num_neighbours = matrix_mean(num_neighbours, nnz);
	*std_num_neighbours = matrix_std_base(num_neighbours, nnz, *avg_num_neighbours);
	free(num_neighbours);
}


__attribute__((unused))
static
double
csr_cross_row_similarity(int * R_offsets, int * C, long m, __attribute__((unused)) long n, __attribute__((unused)) long nnz, long window_size)
{
	int num_threads = omp_get_max_threads();
	double t_row_similarity[num_threads];
	double total_row_similarity;
	long total_num_non_empty_rows = 0;
	#pragma omp parallel
	{
		int tnum = omp_get_thread_num();
		long i, j, k, k_s, k_e, l;
		long degree, num_similarities, column_diff;
		double row_similarity = 0;
		long num_non_empty_rows = 0;
		#pragma omp for schedule(static)
		for (i=0;i<m;i++)
		{
			degree = R_offsets[i+1] - R_offsets[i];
			if (degree <= 0)
				continue;
			for (l=i+1;l<m;l++)       // Find next non-empty row.
				if (R_offsets[l+1] - R_offsets[l] > 0)
					break;
			num_non_empty_rows++;
			if (l < m)
			{
				k_s = R_offsets[l];
				k_e = R_offsets[l+1];
				num_similarities = 0;
				k = k_s;
				for (j=R_offsets[i];j<R_offsets[i+1];j++)
				{
					while (k < k_e)
					{
						column_diff = C[k] - C[j];
						if (labs(column_diff) <= window_size)
						{
							num_similarities++;
							break;
						}
						if (column_diff <= 0)
							k++;
						else
							break;   // went outside of area to examine
					}
				}
				row_similarity += ((double) num_similarities) / degree;
			}
		}
		__atomic_store(&t_row_similarity[tnum], &row_similarity, __ATOMIC_RELAXED);
		__atomic_fetch_add(&total_num_non_empty_rows, num_non_empty_rows, __ATOMIC_RELAXED);
	}
	if (total_num_non_empty_rows == 0)
		return 0;
	total_row_similarity = 0;
	for (long i=0;i<num_threads;i++)
		total_row_similarity += t_row_similarity[i];
	return total_row_similarity / total_num_non_empty_rows;
}


__attribute__((unused))
static
void
csr_cross_row_similarity2(int * R_offsets, int * C, long m, __attribute__((unused)) long n, __attribute__((unused)) long nnz, long window_size,
	double * avg_cross_row_similarity, double * std_cross_row_similarity, double * min_cross_row_similarity, double * max_cross_row_similarity)
{
	double * cross_row_similarity_tmp = (typeof(cross_row_similarity_tmp)) malloc(m * sizeof(*cross_row_similarity_tmp));
	long total_num_non_empty_rows = 0;
	#pragma omp parallel
	{
		long i, j, k, k_s, k_e, l;
		long degree, num_similarities, column_diff;
		long num_non_empty_rows = 0;
		#pragma omp for schedule(static)
		for (i=0;i<m;i++)
		{
			cross_row_similarity_tmp[i] = 0;

			degree = R_offsets[i+1] - R_offsets[i];
			if (degree <= 0)
				continue;
			for (l=i+1;l<m;l++)       // Find next non-empty row.
				if (R_offsets[l+1] - R_offsets[l] > 0)
					break;
			num_non_empty_rows++;
			if (l < m)
			{
				k_s = R_offsets[l];
				k_e = R_offsets[l+1];
				num_similarities = 0;
				k = k_s;
				for (j=R_offsets[i];j<R_offsets[i+1];j++)
				{
					while (k < k_e)
					{
						column_diff = C[k] - C[j];
						if (labs(column_diff) <= window_size)
						{
							num_similarities++;
							break;
						}
						if (column_diff <= 0)
							k++;
						else
							break;   // went outside of area to examine
					}
				}
				cross_row_similarity_tmp[i] = ((double) num_similarities) / degree;
			}

		}
		__atomic_fetch_add(&total_num_non_empty_rows, num_non_empty_rows, __ATOMIC_RELAXED);
	}
	if (total_num_non_empty_rows == 0){
		*avg_cross_row_similarity = 0;
		*std_cross_row_similarity = 0;
		*min_cross_row_similarity = 0;
		*max_cross_row_similarity = 0;
	}
	else{
		// need to keep only those that are !=0 and use right number of "total_num_non_empty_rows", otherwise error compared to initial calculation
		double * cross_row_similarity = (typeof(cross_row_similarity)) malloc(total_num_non_empty_rows * sizeof(*cross_row_similarity));
		long cnt=0;
		for(long i=0; i<m; i++){
			if(cross_row_similarity_tmp[i]>0){
				cross_row_similarity[cnt] = cross_row_similarity_tmp[i];
				cnt++;
			}
		}
		matrix_min_max(cross_row_similarity, total_num_non_empty_rows, min_cross_row_similarity, max_cross_row_similarity);
		*avg_cross_row_similarity = matrix_mean(cross_row_similarity, total_num_non_empty_rows);
		*std_cross_row_similarity = matrix_std_base(cross_row_similarity, total_num_non_empty_rows, *avg_cross_row_similarity);
		free(cross_row_similarity);
	}
	free(cross_row_similarity_tmp);
}


__attribute__((unused))
static
void
csr_ngroups_dis(int * R_offsets, int * C, long m, double * ngroups,
	double * avg_ngroups_size, double * std_ngroups_size, double * min_ngroups_size, double * max_ngroups_size, 
	double * avg_dis, double * std_dis, double * min_dis, double * max_dis)
{	
	long ngroups_total = 0;
	#pragma omp parallel
	{
		long i;
		#pragma omp for schedule(static) reduction(+:ngroups_total)
		for (i=0;i<m;i++)
			ngroups_total+=ngroups[i];
	}

	double * ngroups_size = (typeof(ngroups)) calloc((long)ngroups_total, sizeof(*ngroups));
	double * dis = (typeof(dis)) malloc((long)(ngroups_total-1) * sizeof(*dis ));

	long i, j, cnt=0, last_ci;
	for(i=0; i<m; i++){
		// if(i<17)
		// 	printf("--- row = %ld\n", i);
		long degree = R_offsets[i+1] - R_offsets[i];
		if(degree>0){
			// change of line and we have to consider distance between last element of row i-1 and first element of row i

			// first element, keep its distance from last element of previous row
			if(i>0){
				last_ci = C[R_offsets[i]-1];
				dis[cnt] = fabs((double) (C[R_offsets[i]] - last_ci));
				ngroups_size[cnt]++;
				cnt++;
			}
			// then, proceed with remaining elements of this row
			long prev_ci = C[R_offsets[i]]; 
			for(j=R_offsets[i]+1; j<R_offsets[i+1]; j++){
				ngroups_size[cnt]++;
				if(C[j] > prev_ci+1) { // a new ngroup from now on
					dis[cnt] = fabs((double) (C[j] - prev_ci));
					cnt++;
				}
				prev_ci = C[j];
			}
		}
	}

	// for(i=0;i<22;i++) printf("ngroups_size[%ld] = %ld\tdis[%ld] = %ld\n", i, (long)ngroups_size[i], i, (long)dis[i]);

	matrix_min_max(ngroups_size, ngroups_total, min_ngroups_size, max_ngroups_size);
	*avg_ngroups_size = matrix_mean(ngroups_size, ngroups_total);
	*std_ngroups_size = matrix_std_base(ngroups_size, ngroups_total, *avg_ngroups_size);

	matrix_min_max(dis, ngroups_total-1, min_dis, max_dis);
	*avg_dis = matrix_mean(dis, ngroups_total-1);
	*std_dis = matrix_std_base(dis, ngroups_total-1, *avg_dis);

	free(ngroups_size);
	free(dis);
}

#endif /* MATRIX_UTIL_H */

