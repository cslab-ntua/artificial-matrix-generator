#ifndef RANDOM_H
#define RANDOM_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <complex.h>
#include <math.h>


#define RANDOM_MAX  0x7FFFFFFFL


/*
 * Note:
 *     The functions rand() and srand() conform to SVr4, 4.3BSD, C89, C99, POSIX.1-2001.  The function rand_r() is from POSIX.1-2001.
 *     POSIX.1-2008 marks rand_r() as obsolete.
 *     On older rand() implementations, and on current implementations on different systems, the lower-order bits are much less random than the higher-order bits.
 *     Do not use this function in applications intended to be portable when good randomness is needed.  (Use random(3) instead.)
 *
 * long random(void);
 * void srandom(unsigned seed);
 * char *initstate(unsigned seed, char *state, size_t n);
 * char *setstate(char *state);
 *
 *     The  random() function uses a nonlinear additive feedback random number generator employing a default table of size 31 long integers to return successive pseudo-random numbers in the range from 0 to 2^31 - 1.
 *     The period of this random number generator is very large, approximately 16 * ((2^31) - 1).
 *
 *     The srandom() function sets its argument as the seed for a new sequence of pseudo-random integers to be returned by random().
 *     These sequences are repeatable by calling srandom() with the same seed value.
 *     If no seed value is provided, the random() function is  auto‐ matically seeded with a value of 1.
 *     
 *     The initstate() function allows a state array state to be initialized for use by random().
 *     The size of the state array n is used by initstate() to decide how sophisticated a random number generator it should use—the larger the state array, the better the random numbers will be.
 *     Current "optimal" values for the size of the state array n are 8, 32, 64, 128, and 256 bytes; other amounts will be rounded down to the nearest known amount.
 *     Using less than 8 bytes results in an error.
 *     seed is the seed for the initialization, which specifies a starting point for the random number sequence, and provides for restarting at the same point.
 *     
 *     The  setstate() function changes the state array used by the random() function.
 *     The state array state is used for random number generation until the next call to initstate() or setstate().
 *     state must first have been initialized using initstate() or be the result of a previous call of setstate().
 *
 *     RETURN VALUE
 *         The random() function returns a value between 0 and (2^31) - 1.  The srandom() function returns no value.
 *         The initstate() function returns a pointer to the previous state array.  On error, errno is set to indicate the cause.
 *         On success, setstate() returns a pointer to the previous state array.  On error, it returns NULL, with errno set to indicate the cause of the error.
 *
 * int random_r(struct random_data *buf, int32_t *result);
 * int srandom_r(unsigned int seed, struct random_data *buf);
 * int initstate_r(unsigned int seed, char *statebuf, size_t statelen, struct random_data *buf);
 * int setstate_r(char *statebuf, struct random_data *buf);
 *
 *     The  random_r()  function  is like random(3), except that instead of using state information maintained in a global variable, it uses the state information in the argument pointed to by buf, which must have been previously initialized by initstate_r().  The generated
 *     random number is returned in the argument result.
 *
 *     The srandom_r() function is like srandom(3), except that it initializes the seed for the random number generator whose state is maintained in the object pointed to by buf, which must have been previously initialized by initstate_r(), instead of  the  seed  associated
 *     with the global state variable.
 *
 *     The  initstate_r()  function  is  like  initstate(3)  except  that it initializes the state in the object pointed to by buf, rather than initializing the global state variable.
 *     Before calling this function, the buf.state field must be initialized to NULL.
 *     The initstate_r() function records a pointer to the statebuf argument inside the structure pointed to by buf.
 *     Thus, statebuf should not be deallocated so long as buf is still in use.
 *     (So, statebuf should typically be allocated as a static variable, or allocated on the heap using malloc(3) or similar.)
 *
 *     The  setstate_r()  function is like setstate(3) except that it modifies the state in the object pointed to by buf, rather than modifying the global state variable.  state must first have been initialized using initstate_r() or be the result of a previous call of set‐
 *     state_r().
 *
 *     RETURN VALUE
 *         All of these functions return 0 on success.  On error, -1 is returned, with errno set to indicate the cause of the error.
 *
 * struct random_data
 * {
 * 	int32_t *fptr;               // Front pointer.
 * 	int32_t *rptr;               // Rear pointer.
 * 	int32_t *state;              // Array of state values.
 * 	int rand_type;               // Type of random number generator.
 * 	int rand_deg;                // Degree of random number generator.
 * 	int rand_sep;                // Distance between front and rear.
 * 	int32_t *end_ptr;            // Pointer behind state table.
 * };
 *
 *
 * RANDOM_MAX     = 00000000 00000000 00000000 00000000 01111111 11111111 11111111 11111111
 * RANDOM_MAX + 1 = 00000000 00000000 00000000 00000000 10000000 00000000 00000000 00000000
 *
 * When calling 'random()' / 'random_r()' consecutively we need to take into
 * account instruction reordering if we want reproducible results.
 * E.g.:
 *     if
 *          next(random()) -> 1
 *          next(random()) -> 2
 *     then
 *         x = random() + random() * I
 *     could be one of the two:
 *         1 + 2 * I
 *     or 
 *         2 + 1 * I
 *     because no barriers are inserted.
 *
 * This might also happen if we use '#pragma GCC ivdep' for the initialization loop of an array.
 */


struct Random_State {
	size_t statelen;
	char * statebuf;
	struct random_data * buf;
	unsigned int seed;
};


static inline
struct Random_State *
random_new(unsigned int seed)
{
	struct Random_State * rs;
	size_t statelen;
	char * statebuf;
	struct random_data * buf;

	statelen = 256;

	buf = (typeof(buf)) malloc(sizeof(*buf));
	memset(buf, 0, sizeof(*buf));

	statebuf = (typeof(statebuf)) malloc(statelen);
	memset(statebuf, 0, statelen);

	initstate_r(seed, statebuf, statelen, buf);         // Before calling this function, the buf.state field must be initialized to NULL.
	srandom_r(seed, buf);

	rs = (typeof(rs)) malloc(sizeof(*rs));
	rs->statelen = statelen;
	rs->statebuf = statebuf;
	rs->buf = buf;
	rs->seed = seed;

	return rs;
}


static inline
void
random_destroy(struct Random_State * rs)
{
	free(rs->statebuf);
	free(rs->buf);
	free(rs);
}


static inline
void
random_reseed(struct Random_State * rs, unsigned int seed)
{
	srandom_r(seed, rs->buf);
}


// Return random longs from a uniform distribution over [min, max).
static inline
long
random_uniform_integer(struct Random_State * rs, long min, long max)
{
	long range = max - min;
	int32_t val;
	long res;
	random_r(rs->buf, &val);
	res = val % range + min;
	return res;
}


// Return random doubles from a uniform distribution over [min, max).
static inline
double
random_uniform(struct Random_State * rs, double min, double max)
{
	int32_t val;
	double res;
	random_r(rs->buf, &val);
	res = ((double) val) / ((double) (RANDOM_MAX + 1L));
	res = (max - min) * res + min;
	return res;
}


// Box-Muller transform (wikipedia)
static inline
double
random_normal(struct Random_State * rs, double mean, double std)
{
	const double two_pi = 2.0 * M_PI;
	double u1, u2;

	// Create two uniform random numbers, make sure u1 is greater than epsilon.
	do
	{
		u1 = random_uniform(rs, 0, 1);
		u2 = random_uniform(rs, 0, 1);
	}
	while (u1 <= DBL_EPSILON);

	// Compute z0 and z1.
	double mag = std * sqrt(-2.0 * log(u1));
	double z0  = mag * cos(two_pi * u2) + mean;
	// double z1  = mag * sin(two_pi * u2) + mean;

	return z0;
}


#endif /* RANDOM_H */
