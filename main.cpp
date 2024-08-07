#include "main.h"

void init() {
    fill();      // заполнение полей phi и pi числами
    init_();     // заполнение матриц ковариаций числами
    init_for_generating_initial_conditions();
}

int main()
{
    // Plan
    // Сгенерировать 2 вектора: phi и pi, компоненты которых распределены по Гауссу, то есть 
    // вектор phi_1, phi_2, ..., phi_N, вектор, распределнный по N-мерному Гауссу
    // вектор pi_1, pi_2, ..., pi_N, то же самое, N-мерный Гаусс

    // Это будут начальные условия для 1+1D уравнения Клейна-Гордона (нелинейного),
    // которое нужно с этими начальными условиями решить

    // С помощью полученного решения можно посчитать нужные нам величины (тензор энергии-импульса, 
    // корреляционные функции и т.д.)

    // Далее методом Монте-Карло необходимо взять интеграл по начальным условиям и посчитать, таким 
    // образом эволюцию во времени выбранной наблюдаемой величины

    //solve_with_cond();

    init();

    average();

    //analytical_solution();

    /*
    FILE* phi_sol;
    fopen_s(&phi_sol, "phi.txt", "w+");
    
    if (phi_sol != 0) {
        for (int i = 0; i < SIZE_X; ++i) {
            fprintf(phi_sol, "%4.3lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf\n", x[i], Corr_function_t_0[i], Corr_function_t_1[i], Corr_function_t_2[i], Corr_function_t_3[i], Corr_function_t_4[i], Corr_function_t_5[i]);
        }
        fclose(phi_sol);
    }
    else {
        printf("phi_sol = 0");
    }

    FILE* phi_sol2;
    fopen_s(&phi_sol2, "phi2.txt", "w+");

    if (phi_sol2 != 0) {
        for (int i = 0; i < SIZE_X; ++i) {
            fprintf(phi_sol2, "%4.3lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf\n", x[i], Corr_function_t_0_0[i], Corr_function_t_1_1[i], Corr_function_t_2_2[i], Corr_function_t_3_3[i], Corr_function_t_4_4[i], Corr_function_t_5_5[i]);
        }
        fclose(phi_sol2);
    }
    else {
        printf("phi_sol2 = 0");
    }

    FILE* phi_sol3;
    fopen_s(&phi_sol3, "phi3.txt", "w+");

    if (phi_sol3 != 0) {
        for (int i = 0; i < SIZE_X; ++i) {
            fprintf(phi_sol3, "%4.3lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf\n", x[i], Corr_function_t_0_0_0[i], Corr_function_t_1_1_1[i], Corr_function_t_2_2_2[i], Corr_function_t_3_3_3[i], Corr_function_t_4_4_4[i], Corr_function_t_5_5_5[i]);
        }
        fclose(phi_sol3);
    }
    else {
        printf("phi_sol3 = 0");
    }*/
    
    /*
    FILE* phi_exact_;
    fopen_s(&phi_exact_, "phi_exact.txt", "w+");

    if (phi_exact_ != 0) {
        for (int i = 0; i < SIZE_X; ++i) {            
            fprintf(phi_exact_, "%4.3lf", x[i]);
            for (int k = 0; k < 11; k++) {
                fprintf(phi_exact_, "%10.3lf", phi_exact[i][k]);
            }
            fprintf(phi_exact_, "\n");
        }  
        fclose(phi_exact_);
    }
    else {
        printf("phi_exact_ = 0");
    }*/
    /*
    FILE* phi_asymt_;
    fopen_s(&phi_asymt_, "phi_asymt.txt", "w+");

    if (phi_asymt_ != 0) {
        for (int i = 0; i < SIZE_X; ++i) {
            fprintf(phi_asymt_, "%4.3lf", x[i]);
            for (int k = 0; k < 11; k++) {
                fprintf(phi_asymt_, "%10.3lf", phi_asymt[i][k]);
            }
            fprintf(phi_asymt_, "\n");
        }
        fclose(phi_asymt_);
    }
    else {
        printf("phi_asymt_) = 0");
    }
    */
    
    FILE* phi_aver_;
    fopen_s(&phi_aver_, "phi_aver.txt", "w+");

    if (phi_aver_ != 0) {
        for (int i = 0; i < SIZE_X; ++i) {
            fprintf(phi_aver_, "%4.3lf", x[i]);
            for (int k = 0; k < 11; k++) {
                fprintf(phi_aver_, "%20.3lf", phi_aver[i][k]);
            }
            fprintf(phi_aver_, "\n");
        }
        fclose(phi_aver_);
    }
    else {
        printf("phi_aver_ = 0");
    }
    
    FILE* energy;
    fopen_s(&energy, "energy.txt", "w+");

    if (energy != 0) {
        for (int i = 0; i < SIZE_T; ++i) {
            fprintf(energy, "%4.3lf%20.3lf\n", time_[i], energy_aver[i]);
        }
        fclose(energy);
    }
    else {
        printf("energy = 0");
    }
    
    unsigned int end_time = clock();
    std::cout << "Time needed: " << end_time / 1000. << " sec" << std::endl;
    
}

