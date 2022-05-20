#ifndef SYST_H
#define SYST_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "rand.h"
#include "bfunc.h"
#include "poly.h"
#include "syst.h"
#include <string.h>


/* Type for a system of polynomials. */
typedef struct syst_s syst_t;

// Create a new equation system of zeros
syst_t *syst_new_zero(int m, int n, int d);

/* Create a random system of m equations of degree at most d in d variables. */
syst_t *syst_new_random(int m, int n, int d);


// Sets a coefficient in a polynomial to a specific value
void syst_set_coeff(syst_t *syst, int j, bvar_t i, bool new_val);

// Free the system
void syst_free(syst_t *syst);

// Return n.
int syst_n(syst_t *syst);

// Return m.
int syst_m(syst_t *syst);

// Return d
int syst_d(syst_t *syst);

/* Returns the prefixed solution of the system 
(note: is only valid if the system was made by syst_new_random)*/
int syst_prefixed_sol(syst_t *syst);


// Return the 'j'th equation of the system.
poly_t *syst_equ(syst_t *syst, int j);

// Print the system.
void syst_print(syst_t *syst);
void syst_print_map(syst_t *syst);

// Check if x is a solution.
bool syst_is_solution(syst_t *syst, bvar_t x);
bool syst_is_solution_first_k(syst_t *syst, bvar_t x, int k);


// Return a copy of the system
syst_t *syst_copy(syst_t *syst);

// Write the polynomial system to a file
void save_syst_to_file(syst_t* syst, char* filename);


/* Create a system of r random linear combinations of the
 * equations of  the system*/
syst_t *syst_rand_lin_comb(syst_t* E, int r);

/* Return a new system equal to the join of 'qsyst' and 'r' random linear polynomials. */
//qsyst_t *qsyst_join_linear_polynomials(qsyst_t *qsyst, int r);


/* Compute a solution of 'qsyst' by bruteforce.
 * Return True if a solution exists and put it in 'sol',
 * False otherwise. */
//int qsyst_bruteforce(qsyst_t *qsyst, bvar_t *solution, int number_of_solutions);

#endif
