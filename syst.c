
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "rand.h"
#include "bfunc.h"
#include "poly.h"
#include "syst.h"
#include <string.h>


struct syst_s
{
    int m; // Number of polynomials
    int n; // Number of variables
    int d; // Maximum degree of polynomial system
    bvar_t prefixed_sol;
    poly_t **poly; // Polynomials 
};

syst_t *syst_new_zero(int m, int n, int d)
{
    syst_t *syst;
    syst = malloc(sizeof(syst_t));
    syst->m = m;
    syst->n = n;
    syst->d = d;
    syst->prefixed_sol = 0;
    syst->poly = malloc(m * sizeof(poly_t*));
    for (int i = 0; i < m; i++){    syst->poly[i] = poly_new_zero(n, d);    }

    return syst;    
}

/* NOTE: old - does not enforce a pre-fixed solution
syst_t *syst_new_random(int m, int n, int d)
{
    syst_t *syst;
    syst = malloc(sizeof(syst_t));
    syst->m = m;
    syst->n = n;
    syst->d = d;
    syst->poly = malloc(m * sizeof(poly_t*));
    for (int i = 0; i < m; i++){    syst->poly[i] = poly_new_random(n, d);  }

    return syst;    
}
*/

syst_t *syst_new_random(int m, int n, int d)
{
    syst_t *syst;
    syst = malloc(sizeof(syst_t));
    syst->m = m;
    syst->n = n;
    syst->d = d;
    bvar_t prefixed_sol = bvar_new_random(n);
    syst->prefixed_sol = prefixed_sol;
    syst->poly = malloc(m * sizeof(poly_t*));
    for (int i = 0; i < m; i++){
        syst->poly[i] = poly_new_random(n, d, prefixed_sol); 
    }

    return syst;    
}

void syst_set_coeff(syst_t *syst, int j, bvar_t i, bool new_val){
    poly_set_coeff(syst_equ(syst, j), i, new_val);
}

void syst_free(syst_t *syst)
{
    for (int i = 0; i < syst->m; i++){  poly_free(syst->poly[i]);   }
    free(syst->poly);
    free(syst);
}

int syst_n(syst_t *syst){   return syst->n;}
int syst_m(syst_t *syst){   return syst->m;}
int syst_d(syst_t *syst){   return syst->d;}
int syst_prefixed_sol(syst_t *syst){   return syst->prefixed_sol;}


poly_t *syst_equ(syst_t *syst, int j){  return syst->poly[j];   }

void syst_print(syst_t *syst)
{
    for (int i = 0; i < syst->m; i++){
        printf("(%d) : ", i + 1);
        poly_print(syst->poly[i]);
    }
}

void syst_print_map(syst_t *syst)
{
    for (int i = 0; i < syst->m; i++){
        printf("(%d) : ", i + 1);
        poly_print_map(syst->poly[i]);
    }
}

// TODO: testing a solutions should take about m * binom(n, \downarrow d)
// Which is not true about this procedure?  :( 
bool syst_is_solution(syst_t *syst, bvar_t x)
{
    int i = 0;
    while (i < syst->m && poly_eval(syst->poly[i], x) == 0){    i++;    }
    
    return i >= syst->m;
}


bool syst_is_solution_first_k(syst_t *syst, bvar_t x, int k)
{
    if(k>syst->m){k = syst->m;}

    int i = 0;
    while (i < k && poly_eval(syst->poly[i], x) == 0){    i++;    }
    
    return i >= k;
}


syst_t *syst_copy(syst_t *syst)
{
    syst_t *copy;
    copy = malloc(sizeof(syst_t));
    copy->m = syst->m;
    copy->n = syst->n;
    copy->poly = malloc(syst->m * sizeof(poly_t*));
    for (int i = 0; i < syst->m; i++){  copy->poly[i] = poly_copy(syst->poly[i]);   }
    
    return copy;
}


void save_syst_to_file(syst_t* syst, char* filename){
    
    FILE *file = fopen(filename, "w");
    if(file == NULL){
        printf("#Error: cannot open file\n");
        return;
    }

    char str[1000000];
    sprintf(str, "%d\n", syst_n(syst));
    fwrite(str, 1, strlen(str), file);
    sprintf(str, "%d\n", syst_m(syst));
    fwrite(str, 1, strlen(str), file);
    sprintf(str, "%d\n", syst_d(syst));
    fwrite(str, 1, strlen(str), file);
    sprintf(str, "%d\n", (int)pow(2,syst->n));
    fwrite(str, 1, strlen(str), file);

    poly_t* poly;
    bfunc_t* coeffs;
    char tempstr[100];
    for(int j = 0; j < syst_m(syst); j++){
        sprintf(str, "");
        poly = syst_equ(syst, j);
        coeffs = poly->a;
        for(bvar_t x = (bvar_t)pow(2,syst->n) - 1; x > 0; x--){
            sprintf(tempstr, "%d ", bfunc_get(coeffs, x));
            strcat(str, tempstr);
        }
        sprintf(tempstr, "%d\n", bfunc_get(coeffs, 0));
        strcat(str, tempstr);
        fwrite(str, 1, strlen(str), file);
    }


    fclose(file);
}



syst_t *syst_rand_lin_comb(syst_t* E, int l)
{
    syst_t* E_tilde;
    E_tilde = malloc(sizeof(syst_t));
    E_tilde->m = l;
    E_tilde->n = E->n;
    E_tilde->d = (E->d)*l;
    E_tilde->poly = malloc(l * sizeof(poly_t*));
    
    for (int i = 0; i < l; i++)
    {
        E_tilde->poly[i] = poly_new_zero(E->n, E_tilde->d);
        
        for (int j = 0; j < E->m; j++)
        {
            if (rand_bool()){   poly_add(E_tilde->poly[i], E->poly[j]); }
        }
    }

    return E_tilde;
}
