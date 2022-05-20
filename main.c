
#include "rand.h"
#include "test_dinur2_ours.h"

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "poly.h"
#include "syst.h"
#include <string.h>

#include "rbfunc.h"
#include "dinur2_ours.h"

int main(int argc, char *argv[])
{
    //rand_init(1234321);
    int seed = time(NULL);
    rand_init(seed);


    int n,m,d;

    // for(int i = 10; i <= 18; i++){
    //     n = i, m = (int) n/2, d = 8;
    //     syst_t* E = syst_new_random(m, n, d);

    //     bvar_t prefixed_sol = syst_prefixed_sol(E);

    //     // printf("prefixed solution is valid: %d\n", syst_is_solution(E, prefixed_sol));
        
    //     char filename[100];
    //     sprintf(filename, "datasets/system_%d_%d_%d.txt", n,m,d); 
    //     save_syst_to_file(E, filename);
    // }



    // for(int i = 0; i<10; i++){

    
    //     n = 20, m = 20, d = 2;
    //     syst_t* E = syst_new_random(m, n, d);
    //     int k = round(n/round(2.7*d));
    //     if(k == 0){k = 1;} 
    //     int w = d*(k+1) - k;

    //     bvar_t fes_sols[1 << syst_n(E)];
    //     int nr_sol = 0;

    //     rbfunc_t **ZV = malloc((k + 1) * sizeof(rbfunc_t*));
    //     ZV[0] = rbfunc_new(n - k, w);
    //     int j;
    //     for (j = 1; j < k + 1; j ++){   ZV[j] = rbfunc_new(n - k, w + 1);   }

    //     compute_uvalues(E, k, w, ZV);

    //     for (j = 0; j < k + 1; j++){
    //         rbfunc_free(ZV[j]);
    //     }
    //     free(ZV);
    //     syst_free(E);
    // }


    test_dinur2();    

    

    return 0;
}
