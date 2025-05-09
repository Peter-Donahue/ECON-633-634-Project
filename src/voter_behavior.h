#ifndef VOTER_BEHAVIOR_H
#define VOTER_BEHAVIOR_H

/**
 * Ancillary struct containing parameters needed to calculate posterior of
 * voters. Model parameters are:
 * 	- mu (shock mean),
 * 	- sig (shock standard deviation)
 * 	- nu (upper bound of prior support). 
 * Technical parameters are:
 * 	- beta_l (value at which inner integral d beta_r evaluated at)
 * 	- subint (size of interval partition used in gsl integration routine)
 * 	- h_size (size of history over which integration is performed)
 * 	- integral_type (`l` for mean of beta_l, `r` for mean of beta_r, `m` for
 * 		marginal)
 * 	- gsl_integration_workspace (pointer to gsl integration workspace)
 * Observed histories of policy and economic outcomes stored in gsl vectors:
 * 	- out_obs_h (history of economic outcomes)
 * 	- plc_obs_h (history of policy outcomes)
 */
typedef struct igr_params igr_params;
struct igr_params {
	double beta_l;
	size_t subint;
	gsl_vector* out_obs_h;
	gsl_vector_char* plc_obs_h;
	double mu;
	double sig;
	double nu;
	size_t h_size;
	char integral_type;
	gsl_integration_workspace* giw;
};

/**
 * Functions to perform integration wrt beta_r and beta_l. Integration
 * performed for beta_r first wlog.
 */
double lik_integrand_beta_r(double beta_r, void* igr_p);
double lik_integrand_beta_l(double beta_l, void* igr_p);

/**
 * Wrapper for parameter-wise integration functions. Calculates either
 * posterior means, or marginal.
 */
double clc_integral(char integral_type, void* igr_p);

/**
 * Calculates mean difference between beta_l and beta_r given observed
 * histories of policy and economic outcomes.
 */
double clc_e_diff(igr_params igr_p);

#endif
