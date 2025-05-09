#ifndef STATS_H
#define STATS_H

/* Draws realization of RV u ~ U[0,1]. */
double get_uniform_draw(void);

/**
 * Returns value of pdf of normal distribution centered at `mu`,
 * with standard deviation `sig`, when evaluated at `x`.
 */
double norm_pdf(double x, double mu, double sig);

/* Struct used to store pair of draws from standard normal. */
typedef struct normal_pair normal_pair;
struct normal_pair {
	double z_1;
	double z_2;
};

/**
 * Performs Box-Muller transform to generate a pair of realizations of
 * normal RVs with mean `mu` and standard deviation `sig`.
 */
struct normal_pair get_normal_pair(double mu, double sig);

#endif
