#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <gsl/gsl_matrix.h>

#include "output.h"

const char* mk_fname(char prefix, int K, double sig, bool welf) {

	/* Initialize array to contain file name. */
	static char f_name[25] = {0};

	/* Write parameter-dependent file name to array. */
	if (welf == true) {
		snprintf(f_name, 25, "%c_K%d_sig%1.1f_welf.dat", prefix, K, sig);
	} else {
		snprintf(f_name, 25, "%c_K%d_sig%1.1f.dat", prefix, K, sig);
	}

	return f_name;
}

void write_histories(gsl_matrix* e_diff_mat, gsl_matrix_char* plc_mat,
		     gsl_matrix_uint* plz_mat,
		     int K, double sig, bool welf) {

	/* Initialize prefix container and file name array. */
	char prefix = 'x';
	char f_name[25] = {0};

	/* Write posterior (e_diff). */
	prefix = 'e';
	strcpy(f_name, mk_fname(prefix, K, sig, welf)); 
	int succ_e = 0;
	printf("Writing posterior simulation. ");

	/* Write matrix of histories to file, print result. */
	FILE* ep = fopen(f_name, "w+");
	succ_e = gsl_matrix_fprintf(ep, e_diff_mat, "%.16f");
	fclose(ep);
	if (!succ_e) {
		printf("Completed.\n");
	} else {
		printf("Failed.\n");
	}

	/* Write policy outcomes. */
	prefix = 'c';
	int succ_c;
	strcpy(f_name, mk_fname(prefix, K, sig, welf)); 
	printf("Writing policy simulation. ");

	/* Write matrix of histories to file, print result. */
	FILE* cp = fopen(f_name, "w+");
	succ_c = gsl_matrix_char_fprintf(cp, plc_mat, "%c");
	fclose(cp);
	if (!succ_c) {
		printf("Completed.\n");
	} else {
		printf("Failed.\n");
	}

	/* Write polarization dummy. */
	prefix = 'z';
	int succ_z = 0;
	strcpy(f_name, mk_fname(prefix, K, sig, welf)); 
	printf("Writing polarization dummy simulation. ");

	/* Write matrix of histories to file, print result. */
	FILE* zp = fopen(f_name, "w+");
	succ_z = gsl_matrix_uint_fprintf(zp, plz_mat, "%d");
	fclose(zp);
	if (!succ_z) {
		printf("Completed.\n");
	} else {
		printf("Failed.\n");
	}
}
