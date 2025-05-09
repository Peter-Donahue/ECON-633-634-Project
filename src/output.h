#ifndef OUTPUT_H
#define OUTPUT_H

/**
 * Creates file name to store simulation data, as char array, based
 * on parameter space configuration used for the simulation.
 */
const char* mk_fname(char prefix, int K, double sig, bool welf);

/**
 * Write simulated history to .dat file, named using parameter space
 * configuration used for the simulation.
 */
void write_histories(gsl_matrix* e_diff_mat, gsl_matrix_char* plc_mat,
		     gsl_matrix_uint* plz_mat,
		     int K, double sig, bool welf);
#endif
