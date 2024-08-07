#pragma once


#include <iostream>
#include <stdio.h>
#include <ctime>
#include "PARAMETRS.h"
extern double x[SIZE_X];
extern double time_[SIZE_T];

extern double phi[SIZE_T][SIZE_X];

extern double Corr_function_t_0[SIZE_X];
extern double Corr_function_t_1[SIZE_X];
extern double Corr_function_t_2[SIZE_X];
extern double Corr_function_t_3[SIZE_X];
extern double Corr_function_t_4[SIZE_X];
extern double Corr_function_t_5[SIZE_X];

extern double Corr_function_t_0_0[SIZE_X];
extern double Corr_function_t_1_1[SIZE_X];
extern double Corr_function_t_2_2[SIZE_X];
extern double Corr_function_t_3_3[SIZE_X];
extern double Corr_function_t_4_4[SIZE_X];
extern double Corr_function_t_5_5[SIZE_X];

extern double Corr_function_t_0_0_0[SIZE_X];
extern double Corr_function_t_1_1_1[SIZE_X];
extern double Corr_function_t_2_2_2[SIZE_X];
extern double Corr_function_t_3_3_3[SIZE_X];
extern double Corr_function_t_4_4_4[SIZE_X];
extern double Corr_function_t_5_5_5[SIZE_X];

extern double phi_aver[SIZE_X][11];
extern double energy_aver[SIZE_T];

// Среднее значение phi, аналитическое решение
extern double phi_exact[SIZE_X][11];

// Асимптотики на малых временах
extern double phi_asymt[SIZE_X][11];

void fill();
void init_();
void init_for_generating_initial_conditions();

void solve_with_cond();

void average();

void analytical_solution();
