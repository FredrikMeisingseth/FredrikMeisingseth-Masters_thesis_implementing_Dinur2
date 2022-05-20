#ifndef BVAR_H
#define BVAR_H

#include <stdbool.h>
#include <stdint.h>

#include <stdio.h>
#include "binomials.h"
#include "bvar.h"
#include "rand.h"
#include "string.h"


/* Type for a vector of Boolean variables.
 * For example, 
 * 'bvar_t x = b1011' means 'x_1 = 1, x_2 = 1, x_3 = 0, x_4 = 1. */
typedef uint64_t bvar_t;

/* Maximum number of Boolean variables that can be stored
 * (and work with!) in bvar_t. */
const static int BVAR_MAX = 8 * sizeof(bvar_t) - 1;

// Returns the pointer of a uniformly random boolean variable of length n
// NOTE: due to the limited size of bvar_t, n must not be over 64
bvar_t bvar_new_random(int n);

/* Print the 'n' Boolean variables 'x'.
 * (For testing/debugging.) */
void bvar_print(bvar_t x, int n);

/* Print the 'n' Boolean variables 'x' and their associated value 'y'.
 * (For testing/debugging.) */
void bvar_print_map(bvar_t x, int n, bool y);


/* A function that assigns to *x the next largest n-bit number that
has the same hamming weight if there is such. If else the current hamming 
weight is not the maximal allowed, then the smallest number with the next
larger hamming weight is assigned instead. The function returns 1 if a new
suitable assignment was found and 0 otherwise.

The function iterates over the whole set W_{max_weight}^{n} in order of increasing
hamming weight and increasing number value. 
*/
bool next_subset(bvar_t *x, int *current_weight, int n, int max_weight);

// Input: x* = GrayCode(i*), and n: the maximum length of x* 
// x* and i* are updated as x* = GrayCode(i* + 1) and i* + 1
// Output: a boolean describing wether or not the new i* is out of range
bool GrayCode(bvar_t* x, int* i, int n);

/* Returns the number of values in a restricted boolean function object of 'n' variables and bound weight w. */
bvar_t bvar_number_of_values(int n, int w);

/* Return sum_{i=0}^w binomial(n, i), that is, the number 
 * of binary words of length 'n' and Hamming weigth <= 'w'. */
bvar_t bvar_sum_binomials(int n, int w);

/*output the the index of b in the set of integers of lenght n and weight <= w*/
bvar_t bvar_get_index(bvar_t b, int n);

/* Returns the Hamming weight of the Boolean vector 'x'. */
int bvar_weight(bvar_t x);

/* Returns first n1 bits of x */
bvar_t bvar_first_bits(bvar_t x, int n1);

/* Returns last n2 bits of x as bvar value */
bvar_t bvar_last_bits(bvar_t x, int n, int n1);

#endif
