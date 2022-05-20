
#include <stdio.h>
#include "binomials.h"
#include "bvar.h"
#include "rand.h"
#include "string.h"



bvar_t bvar_new_random(int n){
    bvar_t x = 0;
    for(int i = 0; i < n; i++){
        x += rand_bool()*(1<<i);
    }
    return x;
}

void bvar_print(bvar_t x, int n){
    char s[BVAR_MAX + 1];
    for (int i = n - 1; i >= 0; i--, x >>= 1){
        s[i] = (x & 1) ? '1' : '0';
    }
    s[n] = '\0';
    
    printf("%s\n", s);    
}

void bvar_print_map(bvar_t x, int n, bool y){
    char s[BVAR_MAX + 1];
    for (int i = n - 1; i >= 0; i--, x >>= 1){
        s[i] = (x & 1) ? '1' : '0';
    }
    s[n] = '\0';
    
    printf("%s -> %d\n", s, y);        
}


// TODO: controlcheck this one
bool next_subset(bvar_t *x, int *current_weight, int n, int max_weight){

    //If the current x is 0 then make it 1 if the maximum hamming weight is larger than 0.
    if (*x == 0){
        if (max_weight > 0){
            *x = 1;
            *current_weight = 1;
            return 1;
        }else{
            return 0;
        }
    }

    bvar_t previous_x = *x;
    
    /* 
    Compute next subset as explained in Item 175 of HAKMEM (https://en.wikipedia.org/wiki/HAKMEM).
    The naming convention we use are taken from a blog post by Kate Rose Morley (https://iamkate.com/code/hakmem-item-175/)
    */
    bvar_t lowest_bit, left_bits, changed_bits, right_bits;
    lowest_bit = (*x) & (-*x);  // NOTE: are we actually allowed to do "-" on an unsigned int?
    left_bits = (*x) + lowest_bit;
    changed_bits = (*x) ^ left_bits;
    right_bits = (changed_bits / lowest_bit) >> 2; 
    *x = (left_bits | right_bits);
    
    if(*x >= (bvar_t)1<<n){*x = 0;}
    /*
    If the new x is smaller than the previous one (meaning that there are no greater x
    with the same hamming weight), then increase the hamming weight if the new weight is not
    more than the maximum weight. If it is then break (return 0). 
    */
    // printf("new x: "); bvar_print(*x,64);
    // printf("old x: "); bvar_print(previous_x, 64);
    // printf("new < old: %d\n", *x < previous_x);
    if (*x < previous_x){    
        (*current_weight)++;
        
        if (*current_weight <= max_weight){
            *x = ((bvar_t)1 << *current_weight) - 1;    // Set x to the smallest int of the new hamming weight
            return 1;
        }else{
            return 0;
        }
    }
    
    return 1;
}


bool GrayCode(bvar_t* x, int* i, int n){
    // ffsl gives the position of the lowest set bit in x
    // The transition function of the Gray function is Gray(i + 1) = Gray(i) ^ e_(b_1(i+1))
    // , which means x = x ^ (1 << b_1(i+1))
    *i += 1; 
    *x = (*x) ^ (1 << (ffsl(*i) - 1));
    return *i < (1 << n); 
}

bvar_t GrayCode_from_bvar(bvar_t i, int n){
    bvar_t res = i ^ (i >> 1);
    if(res >= (1<<n)){ res = bvar_first_bits(res, n); }
    return res; 
}


bvar_t bvar_number_of_values(int n, int w){
    return sum_binomials(n, w);
}

bvar_t bvar_get_index(bvar_t b, int n){
    int b_weight = bvar_weight(b);
    bvar_t x = 0;
    int counter = 0;
    if (b){
        for (int i = 0; i < n; ++i) {
            if ((b >> i) & 1){
                counter += 1;
                x += binomial(i, counter);
            }
        }
        return x + bvar_number_of_values(n, b_weight - 1);
    }else{
        return 0;
    }
}

int bvar_weight(bvar_t x){  return __builtin_popcount(x);   }

// NOTE: first in this context mean the LEAST significant bits
bvar_t bvar_first_bits(bvar_t x, int k){   return x & ((1 << k) - 1); }

// NOTE: first in this context mean the LEAST significant bits
/* Returns last n - k bits of x, where x is an bvar object of lenght n */
bvar_t bvar_last_bits(bvar_t x, int n, int k){
    return (x & ((1 << (n - k)) - 1) << k) >> k;
}
