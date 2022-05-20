
#include <stdlib.h>
#include <stdio.h>
#include "bfunc.h"
#include "rbfunc.h"
#include "dinur2_ours.h"

// NOTE: This is the old (non-fes-like version)
// void compute_uvalues(syst_t *E, int k, int w, rbfunc_t **ZV)
// {
//     int n = syst_n(E);
//     bvar_t x, x2, y, z, x1 = 0;
//     int h = 0;
//     do{
//         for (x2 = 0; x2 < (1 << k); x2++)
//         {
//             x = (x1 << k) | x2;
//             if (syst_is_solution(E, x)){
//                 z = bvar_first_bits(x, k);
//                 y = bvar_last_bits(x, n, k);
//                 if (bvar_weight(y) <= w){
//                     rbfunc_xor(ZV[0], y, 1);
//                 }
                
//                 for (int j = 1; j < k + 1; j++){
//                     if ((z & (1 << (k - j))) == 0){
//                         rbfunc_xor(ZV[j], y, 1);
//                     }
//                 }
//             }
//         }
//     }while (next_subset(&x1, &h, n - k, w + 1));
// }

// NOTE: this procedure is not proper FES, amongst other things is not even 
//       the enumeration FES-like
void adapted_FES(syst_t* E, int k, int w, bvar_t* sols, int* nr_sol, int internal_k){
    
    int n = syst_n(E); int d = syst_d(E); int m = syst_m(E);


    bvar_t x, y = 0, z = 0;
    int i = 0, h = 0, sol_index = 0;
    // int current_sol_array_size = (1<<internal_k);

    // bvar_t* tmp_array;
    do{
        i = 0;
        z = 0;
        do{
            x = (y << k) | z;
            if (syst_is_solution_first_k(E, x, internal_k)){
                sols[sol_index] = x;
                sol_index++;
                // if(sol_index >= current_sol_array_size){
                //     current_sol_array_size += (1<<internal_k);
                //     tmp_array = realloc(sols, current_sol_array_size*sizeof(bvar_t));
                //     if (tmp_array == NULL) {
                //         perror("Error allocating memory.\n");
                //         exit(1);
                //     } else {
                //         sols = tmp_array;
                //     }
                // }
            }
        }while(GrayCode(&z, &i, k));
    }while(next_subset(&y, &h, n - k, w + 1));


    int sol_index_new;
    for(; internal_k < m; internal_k++ ){
        sol_index_new = 0;
        for(i = 0; i < sol_index; i++){
            if(poly_eval(syst_equ(E, internal_k), sols[i]) == 0){
                sols[sol_index_new] = sols[i];
                sol_index_new++;
            }
        }
        sol_index = sol_index_new;

    }

    *nr_sol = sol_index;
}

// This is the new version, that actually does (a bit more) FES-like bruteforce
// NOTE: We do naive enumeration and evaluation but atleast do early-abort as in FES
void compute_uvalues(syst_t *E, int k, int w, rbfunc_t **ZV)
{
    int n = syst_n(E); int d = syst_d(E); int m = syst_m(E);
    int internal_k = d*(int)(log2(n)); 
    if(internal_k > m){
        internal_k = (int)log2(m);      //NOTE: this choice should be better motivated
    }
    
    // Note: The expected needed size for fes_sols is 2^(n - internal_k), which means 
    // that this initial size is much too small for small problem sizes, it has however 
    // been found to be a reasonable tradeoff for the range of sizes we can solve
    // bvar_t* fes_sols = malloc((1 << internal_k) * sizeof(bvar_t));
    bvar_t* fes_sols = malloc((1 << n) * sizeof(bvar_t));

    int nr_sol = 0;
    adapted_FES(E, k, w, fes_sols, &nr_sol, internal_k);

    // for(int ii = 0; ii < nr_sol; ii++){
    //     printf("x:");
    //     bvar_print(fes_sols[ii], n);
    // }

    bvar_t x, y, z;
    
    for(int i = 0; i<nr_sol; i++){
        x = fes_sols[i];
        // printf("x:");
        // bvar_print(x, n);
        z = bvar_first_bits(x, k);
        y = bvar_last_bits(x, n, k);
        if (bvar_weight(y) <= w){
            rbfunc_xor(ZV[0], y, 1);
        }
        
        for (int j = 1; j <= k; j++){
            if ((z & (1 << (k - j))) == 0){
                rbfunc_xor(ZV[j], y, 1);
            }
        }
    }

    free(fes_sols);

}
/**/



void output_potential_solutions(syst_t *E, int k, int w, bfunc_t **current_potential_outputs)
{
    int n = syst_n(E);

    rbfunc_t **ZV = malloc((k + 1) * sizeof(rbfunc_t*));
    ZV[0] = rbfunc_new(n - k, w);
    int j;
    for (j = 1; j < k + 1; j ++){   ZV[j] = rbfunc_new(n - k, w + 1);   }
    
    
    compute_uvalues(E, k, w, ZV);
    
    for (j = 0; j < k + 1; j++){    rbfunc_moebius_transform(ZV[j]);   }
    
    bfunc_t **evals;
    evals = malloc((k + 1) * sizeof(bfunc_t*));

    for (j = 0; j < k + 1; j++){
        evals[j] = rbfunc_to_bfunc(ZV[j]);
        rbfunc_free(ZV[j]);
        bfunc_moebius_transform(evals[j]);
    }
    free(ZV);

    bvar_t y;
    for (y = 0; y < (1 << (n - k)); y++){
        if (bfunc_get(evals[0], y)){
            // printf("y_hat:"); bvar_print(y,n);
            bfunc_set(current_potential_outputs[0], y, 1);
            for (int i = 1; i < k + 1; i++){
                bfunc_set(current_potential_outputs[i], y, !bfunc_get(evals[i], y));
            }
        }
    }

    for (j = 0; j < k + 1; j++){    bfunc_free(evals[j]);   }
    free(evals);

    /*
    */
}


int dinur2_solve(syst_t *E, int k, int dF, bvar_t *solution){

    int n = syst_n(E);    
    int l = k+1;
    int w = dF - k;

    int max_iterations = 200;

    bfunc_t **potential_solutions_list[max_iterations];
    bool c = 0;

    for (int iter = 0; iter < max_iterations; iter++){
        
        bfunc_t **current_potential_outputs;
        current_potential_outputs = malloc((k + 1) * sizeof(bfunc_t*));
        for (int i = 0; i < k + 1; i++){
            current_potential_outputs[i] = bfunc_new(n - k);
        }

        
        syst_t *E_tilde = syst_rand_lin_comb(E, l);
        
        //printf("l: %d\n", l);
        //printf("E_tilde:\n");
        //syst_print(E_tilde);

        
        output_potential_solutions(E_tilde, k, w, current_potential_outputs);
        potential_solutions_list[iter] = current_potential_outputs;

        
        bvar_t y;
        bvar_t y_max = 1 << (n - k);
        for (y = 0; y < y_max; y++){
            // printf("y_hat:"); bvar_print(y,n);
            // printf("CurrPotSols_0: %d\n", bfunc_get(current_potential_outputs[0], y));
            if (bfunc_get(current_potential_outputs[0], y)){
                for (int it = 0; it < iter; it++){
                    bool solution_found_before = 1;
                    if (bfunc_get(potential_solutions_list[it][0], y) != 1){
                        solution_found_before = 0;
                    }
                    for (int it2 = 1; it2 < k + 1; it2++){
                        if (bfunc_get(current_potential_outputs[it2], y) !=
                            bfunc_get(potential_solutions_list[it][it2], y)){
                            solution_found_before = 0;
                        }
                    }


                    if (solution_found_before){
                        bvar_t sol = y << k;
                        for (int it3 = 1; it3 <= k; it3++){
                            sol |= bfunc_get(current_potential_outputs[it3], y) << (k - it3);
                        }

                        // printf("Proposed: ");
                        // bvar_print(sol, n);

                        if (syst_is_solution(E, sol)){
                            //The proposed solution was correct! Assign the solution variable to it.
                            *solution = sol;
                            c = 1;

                            // Cleanup
                            int i1, i2;
                            for (i1 = 0; i1 < iter + 1; i1++){
                                for (i2 = 0; i2 < k + 1; i2++){
                                    bfunc_free(potential_solutions_list[i1][i2]);
                                }
                            }

                            for (i1 = 0; i1 < iter + 1; i1++){
                                free(potential_solutions_list[i1]);
                            }

                            syst_free(E_tilde);
                            
                            // End the algorithm
                            return c;
                        }else{
                            // The proposed solution was not correct, go to next proposed solution.
                            break;
                        }
                    }
                }
            }
        }
        /**/
        syst_free(E_tilde);
    }

    // Cleanup in the case that no solution was found
    int i1, i2;
    for (i1 = 0; i1 < max_iterations; i1++)
    {
        for (i2 = 0; i2 < k + 1; i2++)
        {
            bfunc_free(potential_solutions_list[i1][i2]);
        }
    }

    for (i1 = 0; i1 < max_iterations; i1++)
    {
        free(potential_solutions_list[i1]);
    }
    /**/

    return c;
}











// typedef struct FES_state_s FES_state_t;
// struct FES_state_s{
//     bvar_t x;
//     bool y;
//     int i;
// }


// FES_state_t FES_init(poly_t poly, int k_0, bvar_t x_0, int n){
//     FES_state_t state;
//     state->x = x_0;
//     state->y = poly_eval(x_0);
//     state->i = 0;

//     for(int k = 0; i < (1<< n); i++){ 
//         x_0_prim = x_0 ^ GrayCode_from_bvar(1 << (k + k_0), n);
//         // TODO: write a function for getting to the polynomial derivative!
//         poly_t derivative_poly = poly_derivative(poly, k + k_0);       
//         D[k] = FES_init(derivative_poly, k + k_0 + 1, x_0_prim);
//     }

//     return(state);
// }

// void FES_next(FES_state_t* state){

//     state->i += 1; 
//     int k = (ffsl(state->i)
//     (state->x) = (state->x) ^ (1 << k - 1));
//     state->y ^= D[k]->y
//     FES_next(D[k])
// }


// void FES_zeroes(poly_t poly){
//     FES_state_t state = FES_init(poly, 0, 0);

// }




// NOTE: this procedure is not proper FES, amongst other things is not even 
//       the enumeration FES-like
void FES(syst_t* E, bvar_t* sols, int* nr_sol, int internal_k){
    
    int n = syst_n(E); int d = syst_d(E); int m = syst_m(E);
    // int current_sol_array_size = 1<<internal_k;

    bvar_t x = 0;
    int i = 0, h = 0, sol_index = 0;

    // bvar_t* tmp_array;
    do{
        if (syst_is_solution_first_k(E, x, internal_k)){
            sols[sol_index] = x;
            sol_index++;
            // if(sol_index >= current_sol_array_size){
            //     current_sol_array_size += 1<<internal_k;
            //     tmp_array = realloc(sols, current_sol_array_size*sizeof(bvar_t));
            //     if (tmp_array == NULL) {
            //         perror("Error allocating memory.\n");
            //         exit(1);
            //     } else {
            //         sols = tmp_array;
            //     }
            // }
        }
    }while(GrayCode(&x, &i, n));


    int sol_index_new;
    for(; internal_k < m; internal_k++ ){
        sol_index_new = 0;
        for(i = 0; i < sol_index; i++){
            if(poly_eval(syst_equ(E, internal_k), sols[i]) == 0){
                sols[sol_index_new] = sols[i];
                sol_index_new++;
            }
        }
        sol_index = sol_index_new;
    }

    *nr_sol = sol_index;
}
