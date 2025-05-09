#ifndef SIMUL_H
#define SIMUL_H

/**
 * Struct to store state of polity at any given point in time; plc is
 * implemented policy, plz is polarization dummy.
 */
typedef struct polity_state polity_state;
struct polity_state {
	char plc;
	bool plz;
};

/**
 * Struct to store histories of voters' posterior, policy history, and
 * polarization history.
 */
typedef struct history history;
struct history {
	gsl_vector* e_diff_hist;
	gsl_vector_char* plc_hist;
	gsl_vector_uint* plz_hist;
};

/**
 * Generates economic outcome using realization of shock `eps` and
 * policy that is implemented `plc`, given true value of parameters
 * `beta_l_star` and `beta_r_star`.
 */
double mk_out(double eps, int plc, double beta_l_star, double beta_r_star);

/* Generates polity status, i.e. policy outcome and indicator for
 * polarization, given model parameters and voters' posterior.
 */
struct polity_state mk_plt(double e_diff, double zeta, double alpha);

/**
 * Simulates full `T`-long history, given parameter `K`, which identifies memory
 * length, and `sig`, the standard deviation of white noise shocks to
 * economic outcomes.
 */
struct history mk_hist(int K, double sig, size_t T); 

#endif
