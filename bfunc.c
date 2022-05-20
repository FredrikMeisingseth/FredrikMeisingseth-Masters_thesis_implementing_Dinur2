
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include "bfunc.h"

typedef uint64_t bchunk;
const uint64_t bchunk_size = 8 * sizeof(bchunk); 

struct bfunc_s
{
    int        n; /* Number of variables. */
    bvar_t  pow2; /* Number of values. */
    bchunk *data; /* Truth values or monomial coefficients. */
};

bfunc_t *bfunc_new(int n){
    bfunc_t *bfunc;
    
    bfunc = malloc(sizeof(bfunc_t));
    bfunc->n = n;
    bfunc->pow2 = (bvar_t)1 << n;
    bfunc->data = calloc(bfunc->pow2 / bchunk_size + 1, sizeof(bchunk));
    
    return bfunc;
}

void bfunc_free(bfunc_t *bfunc){
    free(bfunc->data);
    free(bfunc);
}

bfunc_t *bfunc_copy(bfunc_t *bfunc){
    bfunc_t *copy = bfunc_new(bfunc->n);
    for (bvar_t x = 0; x < bfunc->pow2; x++){
        bfunc_set(copy, x, bfunc_get(bfunc, x));
    }
    
    return copy;
}

#define BFUNC_GET(bfunc, x) ((bfunc)->data[(x) / bchunk_size] & ((bchunk)1 << ((x) % bchunk_size)))

bool bfunc_get(bfunc_t *bfunc, bvar_t x){
    return BFUNC_GET(bfunc, x);
}

#define BFUNC_SET(bfunc, x) ((bfunc)->data[(x) / bchunk_size] |= ((bchunk)1 << ((x) % bchunk_size)))
#define BFUNC_CLEAR(bfunc, x) ((bfunc)->data[(x) / bchunk_size] &= ~((bchunk)1 << ((x) % bchunk_size)))

void bfunc_set(bfunc_t *bfunc, bvar_t x, bool y)
{
    if (y){ BFUNC_SET(bfunc, x);    }
    else{   BFUNC_CLEAR(bfunc, x);  }
}

void bfunc_add(bfunc_t *bfunc, bvar_t x, bool y){
    bfunc_set(bfunc, x, bfunc_get(bfunc, x) ^ y);
}

void bfunc_print(bfunc_t *bfunc){
    for (int i = 0; i< bfunc->pow2 / bchunk_size + 1; i++){
        bvar_t x = bfunc->data[i];
        bvar_print(x, bchunk_size);
    }
}

void bfunc_print_map(bfunc_t *bfunc){
    for (int i = 0; i< bfunc->pow2 / bchunk_size; i++){
        bvar_t x = bfunc->data[i];
        bvar_print_map(x, bfunc->n, BFUNC_GET(bfunc, x));
    }
}



void bfunc_and(bfunc_t *bfunc1, bfunc_t *bfunc2) {
    for (int i = 0; i < bfunc1->pow2 / bchunk_size + 1; i++){
        bfunc1->data[i] &= bfunc2->data[i];
    }
}

void bfunc_xor(bfunc_t *bfunc1, bfunc_t *bfunc2) {
    for (int i = 0; i < bfunc1->pow2 / bchunk_size + 1; i++){
        bfunc1->data[i] ^= bfunc2->data[i];
    }
}

#define BFUNC_XOR(bfunc, x, y) ((bfunc)->data[(x) / bchunk_size] ^= ((y) ? ((bchunk)1 << ((x) % bchunk_size)) : 0))

/* NOTE: this is the old (standard, non-recursive) implementation
void bfunc_moebius_transform(bfunc_t *bfunc){
    for (int i = 0; i < bfunc->n; i++){
        for(bvar_t x = 0; x < bfunc->pow2; x++){
            if (x & ((bvar_t)1 << i))            {
                BFUNC_XOR(bfunc, x, BFUNC_GET(bfunc, x ^ ((bvar_t)1 << i)));
            }
        }
    }
}
*/

void bfunc_moebius_transform(bfunc_t *bfunc){
    for (int i = 0; i < bfunc->n; i++){
        for(bvar_t x = 0; x < bfunc->pow2; x++){
            if (x & ((bvar_t)1 << i))            {
                BFUNC_XOR(bfunc, x, BFUNC_GET(bfunc, x ^ ((bvar_t)1 << i)));
            }
        }
    }
}
