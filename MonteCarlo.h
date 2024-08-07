#pragma once

#include "PARAMETRS.h"
#include <stdio.h>
#include <math.h>

extern double x[SIZE_X];

extern double phi[SIZE_T][SIZE_X];

extern double Energy[SIZE_T];

void solve_with_cond();
extern double phi_averaged_with_eta;

double Wigner_func();

extern int time_moments[11];

double F(int, int);
