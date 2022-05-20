#ifndef POLY_H
#define POLY_H

#include "bfunc.h"
#include "bvar.h"

/*
Written by Fredrik Meisingseth as a part of his master's thesis at Chalmers and TU Graz. 
The implementation is an extension of the work by Barbero et al (INSERT REF!!)
*/



// NOTE: is it smarter to define the polynomial by its truth table??
/* 
Type for a generic (Boolean) polynomial of arbitrary (but fixed) degree. 
*/
typedef struct poly_s poly_t;

/* Structure for a generic (Boolean) polynomial of arbitrary (but fixed) degree
 * 
 * P(x) = sum_{u in F_2^n} a_{u} x^u,  a_{i,j}, b in F_2
 * 
 * in the variables x_1, ..., x_n (note that x_i^2 = x_i).
 */
struct poly_s
{
    int n; // Number of variables.
    int d; // The degree of the polynomial
    bfunc_t *a; // TODO: store this in a more efficient structure?
};


// Create a zero polynomial in n variables.
poly_t *poly_new_zero(int n, int d);

// Create a random polynomial in n variables and a degree of at most d.
poly_t *poly_new_random(int n, int d, bvar_t prefixed_sol);

// Set the i:th coefficient of the polynomial to new_val 
void poly_set_coeff(poly_t* poly, bvar_t i, bool new_val);

// Return a copy of the polynomial.
poly_t *poly_copy(poly_t *poly);

// Free the polynomial.
void poly_free(poly_t *poly);

// Two ways of printing the polynomial of different verbosity.
void poly_print(poly_t *poly);
void poly_print_map(poly_t *poly);

// Evaluate the polynomial at 'x'.
bool poly_eval(poly_t *poly, bvar_t);

// Add the two polynomials.
void poly_add(poly_t *poly1, poly_t *poly2);


// Returns the derivative of poly with respect to x_i
// poly_t* poly_derivative(poly_t* poly, int i);

#endif
