
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "rand.h"
#include "poly.h"
#include "bvar.h"
#include "bfunc.h"

#include <math.h>


poly_t *poly_new_zero(int n, int d){
    poly_t *poly;
    poly = malloc(sizeof(poly_t));
    poly->n = n;
    poly->d = d;
    poly->a = bfunc_new(n);  
    return poly;
}


poly_t *poly_new_random(int n, int d, bvar_t prefixed_sol){
    poly_t *poly = poly_new_zero(n,d);    
    bvar_t x = 1;
    int h = 1;

    do{
        bfunc_set(poly->a, x, rand_bool());
    }while(next_subset(&x, &h, n, d));
    
    bfunc_set(poly->a, 0, poly_eval(poly, prefixed_sol));
    return poly;
}


void poly_set_coeff(poly_t* poly, bvar_t i, bool new_val){
    bfunc_set(poly->a, i, new_val);
}


void poly_free(poly_t *poly){
    bfunc_free(poly->a);
    free(poly);
}


poly_t *poly_copy(poly_t *poly){
    poly_t *copy;
    
    copy = poly_new_zero(poly->n, poly->d);
    copy->a = bfunc_copy(poly->a);
    copy->d = poly->d;
    
    return copy;
}


void poly_print(poly_t *poly){      bfunc_print(poly->a);       }
void poly_print_map(poly_t *poly){  bfunc_print_map(poly->a);   }


// TODO: rewrite this according to Reed-Muller decomposition
/*
bool poly_eval(poly_t *poly, bvar_t x){
    bvar_t u;
    bool r = 0;
    int n = poly->n;
    bvar_t pow2 = (bvar_t)1<<n;
    for (u = 0; u < pow2; u++){
        // If the coefficient is 1, check if all variables in the monomial are 1
        if(bfunc_get(poly->a, u)){
            r ^= ((u & x) == u);            // TODO: double check that this is equivalent with that all vars in monomial are 1 
        }
    }
    return r;
}
*/

// NOTE: This is the attempt of doing that
// or at least by only iterating over nontrivial coefficients...
bool poly_eval(poly_t *poly, bvar_t x){
    
    bvar_t u = 0;
    int h = 0;
    bool r = 0;
    int n = poly->n;
    int d = poly->d;

    do{
        // If the coefficient is 1, check if all variables in the monomial are 1
        if(bfunc_get(poly->a, u)){
            r ^= ((u & x) == u);            // TODO: double check that this is equivalent with that all vars in monomial are 1 
        }
    }while(next_subset(&u, &h, n, d));

    return r;
}



void poly_add(poly_t *poly1, poly_t *poly2){
    bfunc_xor(poly1->a, poly2->a);
}


// TODO: finish this
// poly_t* poly_derivative(poly_t* poly, int i){
//     poly_t* derivative = poly_new_zero(poly->n, poly->d);

//     return derivative;
// }


