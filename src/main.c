#include <stdio.h>
#include <stdlib.h>
#include "wrappers.h"

#define N_SIMULATIONS_COMPSTAT 20000
#define N_PERIODS_COMPSTAT 110


int main(void) {

	/* Set up number of simulations and history length, run simulation. */
	size_t N = N_SIMULATIONS_COMPSTAT;
	size_t T = N_PERIODS_COMPSTAT;

	double lambda = 1.0;
	
	main_compstat(N, T, lambda);

	return EXIT_SUCCESS;
}
