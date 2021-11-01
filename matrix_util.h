#ifndef MATRIX_UTIL_H
#define MATRIX_UTIL_H


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
		min_partial = max_partial = a[0];
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


#endif /* MATRIX_UTIL_H */

