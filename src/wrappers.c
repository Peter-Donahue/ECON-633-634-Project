#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_integration.h>

#include "stats.h"
#include "voter_behavior.h"
#include "simul.h"
#include "output.h"

#define SEED 12131989

void main_compstat(size_t N, size_t T) {

	/* Initialize arrays of parameter values for memory length of voters
	 * and standard deviation of shocks to economic outcomes. Simulation
	 * is carried out for any combination to perform comparative statics.
	 */
	int mem_arr[2] = {
		mem_arr[0] = 5,
		mem_arr[1] = 10
	};

	double sig_arr[3] = {
		sig_arr[0] = 0.2,
		sig_arr[1] = 1.2,
		sig_arr[2] = 2.5
	};

	/* This is for compstats. Set welfare flag used in output data
	 * filenames to false.
	 */
	bool welf = false;

	/* Initialize seed for reproducibility. */
	unsigned int seed = SEED;
	srand(seed);

	for (size_t s = 0; s < 3; s++) {
		for (size_t k = 0; k < 2; k++) {

			printf("Performing simulation for K = %d, sig = %.1f\n",
					mem_arr[k], sig_arr[s]);

			/* Initialize matrices to contain simulation results */
			gsl_matrix* e_diff_mat = gsl_matrix_calloc(N, T);
			gsl_matrix_char* plc_mat = gsl_matrix_char_calloc(N, T);
			gsl_matrix_uint* plz_mat = gsl_matrix_uint_calloc(N, T);

			/* Run N simulations for each configuration of the
			 * parameter space.
			 */
			for (size_t n = 0; n < N; n++) {
				printf("\tSimulating history %zu\n", n);
				struct history history;
				history = mk_hist(mem_arr[k], sig_arr[s], T);
				gsl_matrix_set_row(e_diff_mat, n,
						   history.e_diff_hist);
				gsl_matrix_char_set_row(plc_mat, n,
							history.plc_hist);
				gsl_matrix_uint_set_row(plz_mat, n,
							history.plz_hist);

			}

			/* Write the N simulated history for each variable of
			 * interest, namely expected difference between betas,
			 * electoral outcome, and polarization dummy.
			 */
			write_histories(e_diff_mat, plc_mat, plz_mat,
			 	        mem_arr[k], sig_arr[s], welf);
			
			/* Free memory. */
			gsl_matrix_free(e_diff_mat);
			gsl_matrix_char_free(plc_mat);
			gsl_matrix_uint_free(plz_mat);

		}
	}

}

void main_welfare(size_t N, size_t T) {

	/* Initialize arrays of parameter values for memory length of voters
	 * and standard deviation of shocks to economic outcomes. Simulation
	 * is carried out for any combination to perform comparative statics.
	 */
	int mem_arr[7] = {
		mem_arr[0] = 5,
		mem_arr[1] = 10,
		mem_arr[2] = 15,
		mem_arr[3] = 20,
		mem_arr[4] = 25,
		mem_arr[5] = 30,
		mem_arr[6] = 50,
	};

	double sig_arr[1] = {
		//sig_arr[0] = 0.5,
		sig_arr[0] = 1.0
	};

	/* This is for welfare analysis. Set welfare flag used in output data
	 * filenames to false.
	 */
	bool welf = true;

	/* Initialize seed for reproducibility. */
	unsigned int seed = SEED;
	srand(seed);

	for (size_t k = 0; k < 7; k++) {
		for (size_t s = 0; s < 1; s++) {

			printf("Performing simulation for K = %d, sig = %.1f\n",
					mem_arr[k], sig_arr[s]);

			/* Initialize matrices to contain simulation results */
			gsl_matrix* e_diff_mat = gsl_matrix_calloc(N, T);
			gsl_matrix_char* plc_mat = gsl_matrix_char_calloc(N, T);
			gsl_matrix_uint* plz_mat = gsl_matrix_uint_calloc(N, T);

			/* Run N simulations for each configuration of the
			 * parameter space.
			 */
			for (size_t n = 0; n < N; n++) {
				printf("\tSimulating history %zu\n", n);
				struct history history;
				history = mk_hist(mem_arr[k], sig_arr[s], T);
				gsl_matrix_set_row(e_diff_mat, n,
						   history.e_diff_hist);
				gsl_matrix_char_set_row(plc_mat, n,
							history.plc_hist);
				gsl_matrix_uint_set_row(plz_mat, n,
							history.plz_hist);

			}

			/* Write the N simulated history for each variable of
			 * interest, namely expected difference between betas,
			 * electoral outcome, and polarization dummy.
			 */
			write_histories(e_diff_mat, plc_mat, plz_mat,
					mem_arr[k], sig_arr[s], welf);

			/* Free memory. */
			gsl_matrix_free(e_diff_mat);
			gsl_matrix_char_free(plc_mat);
			gsl_matrix_uint_free(plz_mat);

			}
	}

}
