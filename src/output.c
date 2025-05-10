#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <gsl/gsl_matrix.h>

#include "output.h"

const char* mk_fname(char* prefix, int K, double sig, bool welf) {

	/* Initialize array to contain file name. */
	static char f_name[32] = {0};

	/* Write parameter-dependent file name to array. */
	if (welf == true) {
		snprintf(f_name, 32, "%7s_K%d_sig%1.1f_welf.dat", prefix, K, sig);
	} else {
		snprintf(f_name, 32, "%7s_K%d_sig%1.1f.dat", prefix, K, sig);
	}

	return f_name;
}

void write_histories(gsl_matrix* e_diff_mat, gsl_matrix* plc_mat,
		     gsl_matrix* plz_mat,
		     int K, double sig, bool welf) {

	/* Initialize prefix container and file name array. */
	char prefix[8] = {0}; //'xxxxxxx';
	char f_name[32] = {0};

	/* Write posterior (e_diff). */
	strcpy(prefix, "MedBeta");      
	//prefix = 'MedBeta';
	strcpy(f_name, mk_fname(prefix, K, sig, welf)); 
	int succ_e = 0;
	printf("Writing posterior simulation. ");
	printf("File name: %s\n", f_name);

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
	strcpy(prefix, "CsnPlcy");      
	int succ_c;
	strcpy(f_name, mk_fname(prefix, K, sig, welf)); 
	printf("Writing policy simulation. ");
	printf("File name: %s\n", f_name);

	/* Write matrix of histories to file, print result. */
	FILE* cp = fopen(f_name, "w+");
	succ_c = gsl_matrix_fprintf(cp, plc_mat, "%.16f");
	fclose(cp);
	if (!succ_c) {
		printf("Completed.\n");
	} else {
		printf("Failed.\n");
	}

	/* Write polarization measure. */
	strcpy(prefix, "Polariz");
	int succ_z = 0;
	strcpy(f_name, mk_fname(prefix, K, sig, welf)); 
	printf("Writing polarization measure simulation. ");
	printf("File name: %s\n", f_name);

	/* Write matrix of histories to file, print result. */
	FILE* zp = fopen(f_name, "w+");
	succ_z = gsl_matrix_fprintf(zp, plz_mat, "%.16f");
	fclose(zp);
	if (!succ_z) {
		printf("Completed.\n");
	} else {
		printf("Failed.\n");
	}
}
