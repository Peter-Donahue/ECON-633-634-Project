// #include <stddef.h>
// #include <gsl/gsl_integration.h>
// #include <gsl/gsl_vector.h>

// #include "stats.h"
// #include "voter_behavior.h"

// #define ABS_ERROR 0x1P-12
// #define REL_ERROR 0x1P-10
// #define LEFT  -1
// #define RIGHT 1

// double lik_integrand_beta_r(double beta_r, void* igr_p) {
// 	/* Calculates, for given beta_l, integrand with respect to beta_r.
// 	 * Depending on char element of struct of integration parameters,
// 	 * returns product of likelihoods itself (key `m`) to calculate
// 	 * the marginal distribution, product of likelihoods times beta_r
// 	 * (key `r`) or times beta_l (key `l`).
// 	 */

// 	static const double left = LEFT;

// 	/* Extract integration and model parameters from ancillary struct. */
// 	struct igr_params* params = (struct igr_params *) igr_p;
// 	double beta_l = (params -> beta_l);
// 	double mu = (params -> mu);
// 	double sig = (params -> sig);
// 	size_t h_size = (params -> h_size);
// 	//char integral_type = (params -> integral_type);
// 	double integral_type = (params -> integral_type);


// 	/* Initialize histories used to calculate posterior. */
// 	gsl_vector* out_obs_h = (params -> out_obs_h);
// 	gsl_vector* plc_obs_h = (params -> plc_obs_h);

// 	/* Define variables used to store integrand, outcome, and policy. */
// 	double lik_integrand_beta_r = 1;
// 	double x;
// 	double out;
// 	double plc;

// 	/* Iterate through observed economic and policy outcomes to construct
// 	 * likelihood.
// 	 */
// 	for (size_t h = 0; h < h_size; h++) {
// 		out = gsl_vector_get(out_obs_h, h);
// 		plc = gsl_vector_get(plc_obs_h, h);

// 		if (plc <= 0) {
// 			x = out - beta_l;
// 		} else {
// 			x = out - beta_r;
// 		}

// 		lik_integrand_beta_r *= norm_pdf(x, mu, sig);
// 	}
// 	/* Iterate through observed outcomes and implemented platform positions */
// 	for (size_t h = 0; h < h_size; h++) {
// 		double out     = gsl_vector_get(out_obs_h, h);
// 		double p_impl  = gsl_vector_get(plc_obs_h, h);
// 		double x       = out - p_impl;
// 		lik_integrand_beta_r *= norm_pdf(x, mu, sig);
// 	}

// 	/* Modify product above (or not if calculating marginal) to construct
// 	 * appropriate integrand.
// 	 */
// 	switch(integral_type) {
// 		case 'l': lik_integrand_beta_r = lik_integrand_beta_r * beta_l;
// 			  break;
// 		case 'r': lik_integrand_beta_r = lik_integrand_beta_r * beta_r;
// 			  break;
// 		case 'm': break;
// 	}

// 	return lik_integrand_beta_r;
// }

// double lik_integrand_beta_l(double beta_l, void* igr_p){
// 	/* Performs integration wrt beta_r. */

// 	/* Extract integration and model parameters from ancillary struct. */
// 	struct igr_params* params = (struct igr_params *) igr_p;
// 	size_t subint = (params -> subint);
// 	double nu = (params -> nu);
// 	(params -> beta_l) = beta_l;

// 	/* Extract integration workspace from ancillary struct. */
// 	gsl_integration_workspace* giw = (params -> giw);

// 	double res, err;
// 	gsl_function R;

// 	R.function = &lik_integrand_beta_r;
// 	R.params = params;

// 	/* Integrate wrt beta_r. */
// 	gsl_integration_qags(&R, 0, nu, ABS_ERROR, REL_ERROR,
// 			     subint, giw, &res, &err);

// 	return res;
// }

// double clc_integral(char integral_type, void* igr_p) {
// 	/* Performs integration with respect to beta_l and beta_r. */

// 	/* Extract integration and model parameters from ancillary struct. */
// 	struct igr_params* params = (struct igr_params *) igr_p;
// 	params -> integral_type = integral_type;
// 	size_t subint = (params -> subint);
// 	double nu = (params -> nu);

// 	/* Initialize integration workspace to calculate outer integral. */
// 	gsl_integration_workspace* iw = gsl_integration_workspace_alloc(subint);

// 	double res, err;
// 	gsl_function L;

// 	L.function = &lik_integrand_beta_l;
// 	L.params = igr_p;

// 	/* Integrate wrt beta_l. */
// 	gsl_integration_qags(&L, 0, nu, ABS_ERROR, REL_ERROR,
// 			     subint, iw, &res, &err);

// 	/* Free second integration workspace. */
// 	gsl_integration_workspace_free(iw);

// 	return res;
// }

// double clc_e_diff(igr_params igr_p) {
// 	/* Calculates, given limited economic and policy outcome histories
// 	 * stored in struct igr_p, the posterior belief of voters on the
// 	 * difference between beta_l and beta_r.
// 	 */

// 	double e_diff;

// 	/* Calculate marginal distribution. */
// 	double marginal = clc_integral('m', &igr_p);

// 	/* Calculate (unscaled) expectation termwise. */
// 	double e_beta_r = clc_integral('r', &igr_p);
// 	double e_beta_l = clc_integral('l', &igr_p);

// 	/* Calculate difference in expectations. */
// 	e_diff = (e_beta_l - e_beta_r) / marginal;

// 	return e_diff;
// }
