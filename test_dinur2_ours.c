
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "test_dinur2_ours.h"
#include "dinur2_ours.h"
#include "rand.h"
#include "syst.h"
#include <string.h>

void test_dinur2_solver(syst_t* E, int m, int n, int d, int k, int dF, int repeats, bool generate_system)
{
    int c, iteration;
    for (c = 0, iteration = 0; iteration < repeats; iteration++)
    {
        if(generate_system){    E = syst_new_random(m, n, d);   }
        // printf("E:\n");
        // syst_print(E);
        // printf("pre-fixed solution:");
        //  bvar_print(syst_prefixed_sol(E), n);

        printf("System #%d\n", iteration + 1);

        bvar_t sol;
        bool sol_has_been_found = dinur2_solve(E, k, dF, &sol);

        if (sol_has_been_found){
            printf("Solution found: \n");
            bvar_print(sol, n);
        }else{
            printf("No solution found.\n");
        }
        
        
        if(sol_has_been_found){
            if (syst_is_solution(E, sol)){
                printf("---> Right!\n");
                c++;
            }else{
                printf("---> Wrong!\n");
            }
        }
        printf("\n");
        if(generate_system){    syst_free(E);   }
    }
    
    printf("Found correct sol %d out of %d (%.2lf%%) times \n", c, repeats, 100 * c / (double)repeats);
}


syst_t* readFile(syst_t* E, char* filename){
    /*
    Reads the file and parses it into a polynomial system.

    The layout of the input file must be:
    n
    m
    d
    <length of coefficient rows>
    coefficients in rows (in lex ordering!)
    */

    FILE *file = fopen(filename, "r");
    if(file == NULL){
        printf("#Error: cannot open file\n");
        return -1;
    }

    fseek(file, 0L, SEEK_END);
    int file_size = ftell(file);
    rewind(file);
    char *init = (char*) malloc(sizeof(char) * file_size);
    fread(init, sizeof(char), file_size, file);
    fclose(file);
    //printf("%s",init);

    int start, ix, arg_nr;
    int n, m, d, nr_coeffs;
    int totlen = strlen(init);
    for(start = 0, ix = 0, arg_nr = 0; arg_nr < 4; ++ix){
        if(init[ix] == ' ' || init[ix] == '\n'){
            init[ix] = 0;
            if(arg_nr == 0){n = atoi(&init[start]);}
            if(arg_nr == 1){m = atoi(&init[start]);}
            if(arg_nr == 2){d = atoi(&init[start]);}
            if(arg_nr == 3){nr_coeffs = atoi(&init[start]);}
            start = ix + 1;
            ++arg_nr;
        }
    }
    printf("n: %d, m: %d, d:%d\n",n,m,d);
    E = syst_new_zero(m, n, d);

    bvar_t coeff_idx = nr_coeffs - 1;
    int j = 0;

    // printf("\nfile: %s\n\n", &init[start]);
    for( ; ix < totlen, j<m ; ++ix ){
        // printf("ix: %d, start: %d, coeff_idx: %d, j: %d \n", ix, start, coeff_idx, j);
        if(init[ix] == '\n'){
            init[ix] = 0;
            syst_set_coeff(E, j, coeff_idx, atoi(&init[start]));
            start = ix + 1;
            ++j; coeff_idx = nr_coeffs-1;
        }
        if(init[ix] == ' '){
            init[ix] = 0;
            syst_set_coeff(E, j, coeff_idx, atoi(&init[start]));
            start = ix + 1;
            --coeff_idx;
        }
    }
    
    // printf("E:\n");
    // syst_print(E);
    free(init);
    return E;
}

void test_dinur2()
{
    int m, n, d, k, dF, repeats;

    bool generate_system = 0;
    syst_t* E;
    if(!generate_system){        
        char* inputfile_name = "datasets/system_18_9_4.txt";
        E = readFile(E, inputfile_name);
        n = syst_n(E);
        m = syst_m(E);
        d = syst_d(E);
        // printf("E:\n");
        // syst_print(E);
    }else{
        n = 16;
        m = (int)2*n;
        d = 2; 
    }
    
    // printf("n: %d, m: %d, d: %d\n", n, m ,d);
    // printf("E:\n");
    // syst_print(E);
        
    k = round(n/round(2.7*d));      // Note: this one is not necessarily optimal
    if(k == 0){k = 1;}              // Note: Why do i need this? 
    dF = d*(k+1);                   // Note: could this one BE any lower?
    repeats = 1;
    test_dinur2_solver(E, m, n, d, k, dF, repeats, generate_system);

    
    // int internal_k = d*(int)(log2(n)); 
    // if(internal_k > m){
    //     internal_k = (int)log2(m);      //NOTE: this choice should be better motivated
    // }
    // // bvar_t* fes_sols = malloc((1 << (internal_k)) * sizeof(bvar_t));
    // bvar_t* fes_sols = malloc((1 << n) * sizeof(bvar_t));

    // int nr_sol = 0;
    // FES(E, fes_sols, &nr_sol, internal_k);

    // for(int i = 0; i < nr_sol; i++){
    //     printf("sol %d:", i);
    //     bvar_print(fes_sols[i], n);
    // }


    // if(!generate_system){
    //     syst_free(E);
    // }
    
}
