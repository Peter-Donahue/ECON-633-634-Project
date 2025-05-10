#ifndef WRAPPERS_H
#define WRAPPERS_H

/**
 * Performs main simultion for comparative statics purposes. Iterates through
 * array of predefined values of standard deviation of shocks to economic
 * outcomes and memory lengths, simulates N T-length histories.
 */
void main_compstat(size_t N, size_t T, double lambda);

/**
 * Performs main simulation for welfare analysis. Iterates through more values
 * of memory length, with a longer time series duration.
 */
void main_welfare(size_t N, size_t T, double lambda);

#endif
