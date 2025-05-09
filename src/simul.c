#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stddef.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_integration.h>

#include "simul.h"
#include "stats.h"
#include "voter_behavior.h"

#define LEFT  'l'
#define RIGHT 'r'
#define BETA_L_STAR 3.5
#define BETA_R_STAR 2.5

#define SUBINT_PARTITION_CARDINALITY 0x1P+10

#define ALPHA 0.5
#define ZETA 0.5
#define MU 0.0
#define NU 6

double mk_out(double eps, int plc, double beta_l_star, double beta_r_star) {

	/* Initialize policy flags and economic outcome. */
	static char left = LEFT;
	double out = 0;
	
	/* Economic outcome is parameter corresponding to implemented policy
	 * plus normal shock.
	 */
	if (plc == left) {
		out = beta_l_star + eps;
	} else {
		out = beta_r_star + eps;
	}

	return out;
}

struct polity_state mk_plt(double e_diff, double zeta, double alpha) {

	/* Initialize policy flags. */
	static const char left = LEFT;
	static const char right = RIGHT;
	
	/* Initialize electoral threshold. */
	const double rho = (1 / (2 * zeta * (1 + alpha)));
	
	/* Initialize winning probability for `l` in case of indifference. */
	const double p_win_l = 1/2 + zeta * e_diff;

	/* Initialize polity state. */
	struct polity_state plt;
	plt = (polity_state) {'u', true};

	/* Run elections. In first case left wins, with consensus. In second
	 * case right wins, with consensus. In third case polarization occurs,
	 * and winning policy is stochastic, depends on voters' preference
	 * shock phi.
	 */
	if (e_diff > rho) {
		plt.plc = left;
		plt.plz = false;
	} else if (-e_diff > rho) {
		plt.plc = right;
		plt.plz = false;
	} else {
		/* Polarization occurs, i.e. `L` offers `l` and `R` offers `r`.
		 * Left wins with probability `p_win_l`. To determine electoral
		 * outcome, draw u from U[0,1], and assign victory to `L` iff
		 * u <= p_win_l, since F(u) = u.
		 */
		double u = get_uniform_draw();
		if (u <= p_win_l) {
			plt.plc = left;
		} else {
			plt.plc = right;
		}
		plt.plz = true;
	}

	return plt;
}

struct history mk_hist(int K, double sig, size_t T) {
	
	/* Set model parameters. */
	double mu = MU;
	double zeta = ZETA;
	double alpha = ALPHA;
	double nu = NU;
	static const double beta_l_star = BETA_L_STAR;
	static const double beta_r_star = BETA_R_STAR;

	/* Initialize history vectors*/
	gsl_vector* e_diff_hist = gsl_vector_calloc(T);
	gsl_vector_char* plc_hist = gsl_vector_char_calloc(T);
	gsl_vector_uint* plz_hist = gsl_vector_uint_calloc(T);
	gsl_vector* eps_hist = gsl_vector_calloc(T);
	gsl_vector* out_hist = gsl_vector_calloc(T);

	/* Initialize and draw history of shocks. Since working with pairs of
	 * RVs, draw up to second last period.
	 */
	for (size_t t = 0; t < T - 1; t+= 2) {
		struct normal_pair normal_pair;
		normal_pair = get_normal_pair(mu, sig);
		gsl_vector_set(eps_hist, t, normal_pair.z_1);
		gsl_vector_set(eps_hist, t+1, normal_pair.z_2);
	}
	
	/* Initialize observed policy and outcome histories as views on full
	 * policy and outcome histories. Observed histories increase in length
	 * when simulating, until length K is reached. At that point vector
	 * views start shifting and have fixed length K.
	 */
	gsl_vector_view out_obs_h;
	gsl_vector_char_view plc_obs_h;
	size_t h_size = 1;

	/* Initialization of e_diff = 0 at t=0 implies polity status
	 * is determined by coin toss in that period. After that, use
	 * observed histories to calculate posterior.
	 */
	double e_diff = 0.0;
	double out = 0.0;
	double eps = 0.0;
	struct polity_state plt;

	/* Set up struct containing model parameters and integration workspace
	 * to calculate posterior.
	 */
	size_t subint = SUBINT_PARTITION_CARDINALITY;
	gsl_integration_workspace* giw = gsl_integration_workspace_alloc(subint);
	igr_params igr_p = {.subint=subint, .mu=mu, .sig=sig, .nu=nu,
			    .h_size=h_size, .giw=giw};


	for (size_t t = 0; t < T; t++) {

		/* Calculate posterior difference between beta_l and beta_r,
		 * using observed policy and outcome histories.
		 */
		if (t > 0) {
			igr_p.out_obs_h = &out_obs_h.vector;
			igr_p.plc_obs_h = &plc_obs_h.vector;
			e_diff = clc_e_diff(igr_p);
		}

		/* Run elections given median voter's posterior */
		plt = mk_plt(e_diff, zeta, alpha);

		/* Determine economic outcome, given shock and policy implemented
		 * as result of electoral outcome.
		 */
		eps = gsl_vector_get(eps_hist, t);
		out = mk_out(eps, plt.plc, beta_l_star, beta_r_star);

		/* Store outcomes and posterior in vectors. */
		gsl_vector_set(e_diff_hist, t, e_diff);
		gsl_vector_char_set(plc_hist, t, plt.plc);
		gsl_vector_uint_set(plz_hist, t, plt.plz);
		gsl_vector_set(out_hist, t, out);

		if (h_size <= K) {
			/* Outcome and policy histories used for posterior
			 * formation increase in length until length reaches
			 * maximum memory length, i.e. K.
			 */
			out_obs_h = gsl_vector_subvector(out_hist, 0, h_size);
			plc_obs_h = gsl_vector_char_subvector(plc_hist, 0,
							      h_size);
			/* Update integration parameter and increase length of
			 * history
			 */
			igr_p.h_size = h_size;
			h_size++;

		} else {
			/* When maximum length of observed histories is reached,
			 * histories used to form posteriors are views on full
			 * histories, whose starting point is an offset of the
			 * current `t` of length K.
			 */
			out_obs_h = gsl_vector_subvector(out_hist, t - K + 1, K);
			plc_obs_h = gsl_vector_char_subvector(plc_hist, t - K + 1, K);
		}

	}

	/* Store histories. */
	history history = {.e_diff_hist=e_diff_hist,
			   .plc_hist=plc_hist,
			   .plz_hist=plz_hist};

	/* Clean up memory. */
	gsl_vector_free(eps_hist);
	gsl_integration_workspace_free(giw);

	return history;
}
