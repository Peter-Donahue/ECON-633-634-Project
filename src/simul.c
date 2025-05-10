#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stddef.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_multifit.h>

#include "simul.h"
#include "stats.h"
// #include "voter_behavior.h"

#define LEFT  'l'
#define RIGHT 'r'

#define EXTREME_LEFT -1.0
#define EXTREME_RIGHT 1.0

#define BETA_L_STAR 3.5
#define BETA_R_STAR 2.5

#define SUBINT_PARTITION_CARDINALITY 0x1P+10

#define ALPHA 0.5
#define ZETA 0.5
#define MU 0.0
#define NU 6

#define USE_POLYNOMIAL_REG 1


double mk_out(double eps, double plc, double beta_l_star, double beta_r_star) {
	double out = 0;
	// We are going to interpolate for now between the beta_l_star and beta_r_star
	const double ext_left = EXTREME_LEFT;
	const double ext_right = EXTREME_RIGHT;

	double delta_beta = 0;
	double delta_policy_val = 0;
	double intercept = 0;
	double slope = 0;

	#if USE_POLYNOMIAL_REG
//     	A polynomial with a max at around -0.5
		out = 2 + (-2/3) * plc + (-0.6) * plc * plc + eps;
	#else
		delta_policy_val = ext_right - ext_left;
		delta_beta = beta_r_star - beta_l_star;
		slope = delta_beta / delta_policy_val;
		intercept = beta_l_star + (delta_policy_val / 2) * slope;
		out = slope * plc + intercept + eps;
		//printf("out = %f, slope = %f, plc = %f, intercept = %f, eps = %f\n", out, slope, plc, intercept, eps);
	#endif
	return out;
}

/*
 * Solve: min_p [(p - opt_voter_p)^2 + lambda * (p - I)^2]
 * Analytic solution: p* = (opt_voter_p + lambda * I) / (1 + lambda), clamped to [-1,1]
 */
double choose_position(double opt_voter_p, double I, double lambda) {
    double p = (opt_voter_p + lambda * I) / (1.0 + lambda);
    if (p < -1.0) p = -1.0;
    if (p > +1.0) p = +1.0;
    return p;
}

struct polity_state mk_plt(double opt_voter_p, double zeta, double alpha, double lambda) {
	/* party ideals */
	const double ideal_l = -1.0;
	const double ideal_r = +1.0;
	
    /* compute each party's optimal platform p_l, p_r */
    double p_l = choose_position(opt_voter_p, ideal_l, lambda);
    double p_r = choose_position(opt_voter_p, ideal_r, lambda);

	/* median voter belief = opt_voter_p */
    /* determine winner: which platform is numerically closer to median belief */
    double d_l = fabs(p_l - opt_voter_p);
    double d_r = fabs(p_r - opt_voter_p);
	
	bool left_wins;
	if (d_l < d_r) {             
		left_wins = true;
	} else if (d_r < d_l) {         
		left_wins = false;
	} else {                        
		left_wins = (get_uniform_draw() < 0.5);
	}

	struct polity_state plt;
	plt.plc      = left_wins ? p_l : p_r;
	plt.winner_l = left_wins;
	plt.plz      = fabs(p_l - p_r); /* polarization measure */
	return plt;
}

/* compute optimal policy from quadratic coefficients */
static double optimal_policy(double B1, double B2) {
    if (fabs(B2) < 1e-12)            /* nearly flat */
        return (B1 < 0.0 ?  -1.0 : 1.0);
    double p = B1 / (2.0 * B2);
    if (p < -1.0) p = -1.0;
    if (p >  1.0) p =  1.0;
    return p;
}



struct history mk_hist(int K, double sig, size_t T, double lambda) {
	/* ---- model parameters ---- */
	const double mu   = MU;
	const double zeta = ZETA;
	const double alpha= ALPHA;
	const double beta_l_star = BETA_L_STAR;
	const double beta_r_star = BETA_R_STAR;

	/* ---- storage ---- */
	gsl_vector *opt_voter_p_hist   = gsl_vector_calloc(T);      
	gsl_vector *plc_hist = gsl_vector_calloc(T);       /* implemented platform */
	gsl_vector *plz_hist = gsl_vector_calloc(T); /* winner flag */
	gsl_vector *eps_hist = gsl_vector_calloc(T);       /* iid shocks */
	gsl_vector *out_hist = gsl_vector_calloc(T);       /* realized outcomes */

	/* pre‑draw shocks */
	for (size_t t=0; t<T-1; t+=2) {
		struct normal_pair z = get_normal_pair(mu, sig);
		gsl_vector_set(eps_hist, t,   z.z_1);
		gsl_vector_set(eps_hist, t+1, z.z_2);
	}
	/* rolling views */
	gsl_vector_view out_obs_h;
	gsl_vector_view plc_obs_h;
	size_t h_size = 1; /* grows to K */
	for (size_t t=0; t<T; ++t) {
		/* ---------- voters form belief via OLS y = B·p ---------- */
		//double B1 = 0.0; /* default when no past data */
		//double B2 = 0.0; /* default when no past data */
		//double intercept = 0.0;
		double opt_voter_p = 2 * get_uniform_draw() - 1.0;  /* random in [-1,1] */
		if (t>0) {
			if (h_size <= K) {
				out_obs_h = gsl_vector_subvector(out_hist, 0, h_size);
				plc_obs_h = gsl_vector_subvector(plc_hist, 0, h_size);
			} else {
				out_obs_h = gsl_vector_subvector(out_hist, t-K, K);
				plc_obs_h = gsl_vector_subvector(plc_hist, t-K, K);
			}
			
			size_t n = out_obs_h.vector.size;
			#if USE_POLYNOMIAL_REG
				const size_t p = 3;
			#else
				const size_t p = 2;
			#endif
			if (n>=p) {
				//const size_t p = 1;                        /* # regressors */
				gsl_matrix *X = gsl_matrix_alloc(n, p);
				gsl_vector *y = gsl_vector_alloc(n);
				for (size_t i = 0; i < n; ++i) {
					double pval = gsl_vector_get(&plc_obs_h.vector, i);
					double oval = gsl_vector_get(&out_obs_h.vector, i);
					gsl_matrix_set(X, i, 0, 1.0);    /* intercept */
					gsl_matrix_set(X, i, 1, pval);   /* p     */
					#if USE_POLYNOMIAL_REG
						gsl_matrix_set(X, i, 2, pval*pval); /* p^2     */
					#endif
					gsl_vector_set(y,     i, oval);
				}
				
				gsl_vector *coef = gsl_vector_alloc(p);
				gsl_matrix *cov  = gsl_matrix_alloc(p, p);
				double chisq;
				gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(n, p);
				gsl_multifit_linear(X, y, coef, cov, &chisq, work);
				double intercept = gsl_vector_get(coef, 0);
				double B1 = gsl_vector_get(coef, 1);
				#if USE_POLYNOMIAL_REG
                	double B2 = gsl_vector_get(coef, 2);
				#else
					double B2 = 0.0;
				#endif
                opt_voter_p = optimal_policy(B1, B2);
				//printf("slope = %f, intercept = %f\n", M, intercept);
				gsl_multifit_linear_free(work);
				gsl_matrix_free(X);
				gsl_vector_free(y);
				gsl_vector_free(coef);
				gsl_matrix_free(cov);
			}
		}

		struct polity_state plt = mk_plt(opt_voter_p, zeta, alpha, lambda);
		double eps = gsl_vector_get(eps_hist, t);
		double out = mk_out(eps, plt.plc, beta_l_star, beta_r_star);	

		/* store */
		gsl_vector_set(opt_voter_p_hist, t, opt_voter_p);
		gsl_vector_set(plc_hist, t, plt.plc);
		gsl_vector_set(plz_hist, t, plt.plz);
		gsl_vector_set(out_hist, t, out);

		if (h_size < K) h_size++;
	}

	struct history hist = {
		.e_diff_hist = opt_voter_p_hist,    /* keeping original field name */
		.plc_hist    = plc_hist,
		.plz_hist    = plz_hist
	};

	gsl_vector_free(eps_hist);
	gsl_vector_free(out_hist);
	return hist;
}
