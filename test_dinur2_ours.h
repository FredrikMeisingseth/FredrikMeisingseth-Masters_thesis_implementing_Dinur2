#ifndef TEST_DINUR2_H
#define TEST_DINUR2_H

#include "syst.h"


void test_dinur2_solver(syst_t* E, int m, int n, int d, int k, int dF, int repeats, bool generate_system);
void test_dinur2();
syst_t* readFile(syst_t* E, char* filename);

#endif
