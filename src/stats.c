#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "stats.h"

#define M_PI (3.14159265358979323846)
#define BM_TOLERANCE 1E-6

double get_uniform_draw(void) {

	/* Scale random integer by largest random integer, recast as double. */
	return ((double)rand()/(double)RAND_MAX);
}

double norm_pdf(double x, double mu, double sig) {

	/* Calculate constant factor of density. */
	const long double pi = M_PI;
	const long double inv_sqrt_2pi = 1 / sqrt(2 * pi);

	double a = (x - mu) / sig;

	/* Return density. */
	return inv_sqrt_2pi / sig * exp(-0.5 * a * a);
}

struct normal_pair get_normal_pair(double mu, double sig) {

	/* Uniform draw in ln() must be at least as large as value below. */
	static const double tol = BM_TOLERANCE;
	static const double pi = M_PI;

	/* Draw pair of uniform realizations. Keep drawing until argument of
	 * ln() above threshold.
	 */
	double u_1, u_2;
	do {
		u_1 = get_uniform_draw();
	} while (u_1 <= tol);
	u_2 = get_uniform_draw();

	/* Perform BM transform. */
	double R = sig * sqrt(-2.0 * log(u_1));
	double z_1 = R * cos(2 * pi * u_2) + mu;
	double z_2 = R * sin(2 * pi * u_2) + mu;

	/* Initialize struct to store results, store results, return them. */
	normal_pair z_pair;
	z_pair = (normal_pair) {z_1, z_2};

	return z_pair;
}
