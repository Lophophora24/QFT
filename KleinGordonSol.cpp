#include "KleinGordonSol.h"

double x[SIZE_X];
double time_[SIZE_T];

double Energy[SIZE_T];

double phi[SIZE_T][SIZE_X];

int time_moments[11] = {
    (int)(0 * SIZE_T),
    (int)(0.1 * SIZE_T),
    (int)(0.2 * SIZE_T),
    (int)(0.3 * SIZE_T),
    (int)(0.4 * SIZE_T),
    (int)(0.5 * SIZE_T),
    (int)(0.6 * SIZE_T),
    (int)(0.7 * SIZE_T),
    (int)(0.8 * SIZE_T),
    (int)(0.9 * SIZE_T),
    (int)(0.95 * SIZE_T),
};

// Среднее значение phi, аналитическое решение
double phi_exact[SIZE_X][11];

double phi_asymt[SIZE_X][11];

void fill()
{
    for (int i = 0; i < SIZE_X; ++i)
        x[i] = x1 + i * h;
    for (int i = 0; i < SIZE_T; ++i)
        time_[i] = t1 + i * t;
}

void fill_phi()
{
    for (int i = 0; i < SIZE_X; ++i) {
        for (int j = 0; j < SIZE_T; ++j) {
            phi[j][i] = 0;
            Energy[j] = 0;
        }
        phi[0][i] = phi_0[i];
        phi[1][i] = phi[0][i] + pi_0[i] * t;
    }
}

double F(int j, int i)
{
    //return m * m * phi[j][i] * phi[j][i] / 2 + g * phi[j][i] * phi[j][i] * phi[j][i] * phi[j][i] / 4;

    return pow(m, 4) / (6 * g) * cos(sqrt(6 * g) / m * phi[j][i]);

    //return m * m * phi[j][i] * phi[j][i] / 2 + g * phi[j][i] * phi[j][i] * phi[j][i] * phi[j][i] / 4 + g * g / (20 * m * m) * pow(phi[j][i], 6);
}

void Sinh_Gordon()
{
    for (int j = 1; j < SIZE_T - 1; ++j) {
        for (int i = 0; i < SIZE_X; ++i) {
            if (i == 0) {
                phi[j + 1][i] = r * r * phi[j][SIZE_X - 2] + 2 * (1 - r * r) * phi[j][i] + r * r * phi[j][i + 1] - phi[j - 1][i] - t * t * pow(m, 3) / sqrt(6 * g) * sinh(sqrt(6 * g) / m * phi[j][i]);
            }
            else
            if (i == (SIZE_X - 1)) {
                phi[j + 1][i] = r * r * phi[j][i - 1] + 2 * (1 - r * r) * phi[j][i] + r * r * phi[j][1] - phi[j - 1][i] - t * t * pow(m, 3) / sqrt(6 * g) * sinh(sqrt(6 * g) / m * phi[j][i]);
            }
            else
            {
                phi[j + 1][i] = r * r * phi[j][i - 1] + 2 * (1 - r * r) * phi[j][i] + r * r * phi[j][i + 1] - phi[j - 1][i] - t * t * pow(m, 3) / sqrt(6 * g) * sinh(sqrt(6 * g) / m * phi[j][i]);
            }
        }
    }
}

void Strauss_Vazquez()
{
    for (int j = 1; j < SIZE_T - 1; ++j) {
        for (int i = 0; i < SIZE_X; ++i) {
            if (i == 0) {
                double recent = 100;
                phi[j + 1][i] = phi[j][i];
                do {
                    recent = phi[j + 1][i];
                    if (phi[j + 1][i] == phi[j - 1][i]) {
                        phi[j + 1][i] = r * r * phi[j][SIZE_X - 2] + 2 * (1 - r * r) * phi[j][i] + r * r * phi[j][i + 1] - phi[j - 1][i] - t * t * (m * m * phi[j + 1][i] + g * pow(phi[j + 1][i], 3));
                    }
                    else {
                        phi[j + 1][i] = r * r * phi[j][SIZE_X - 2] + 2 * (1 - r * r) * phi[j][i] + r * r * phi[j][i + 1] - phi[j - 1][i] - t * t * (F(j + 1, i) - F(j - 1, i)) / (phi[j + 1][i] - phi[j - 1][i]);
                    }

                } while (abs(recent - phi[j + 1][i]) > 1e-5);
            }
            else
            if (i == (SIZE_X - 1)) {
                double recent = 100;
                phi[j + 1][i] = phi[j][i];
                do {
                    recent = phi[j + 1][i];
                    if (phi[j + 1][i] == phi[j - 1][i]) {
                        phi[j + 1][i] = r * r * phi[j][i - 1] + 2 * (1 - r * r) * phi[j][i] + r * r * phi[j][1] - phi[j - 1][i] - t * t * (m * m * phi[j + 1][i] + g * pow(phi[j + 1][i], 3));
                    }
                    else {
                        phi[j + 1][i] = r * r * phi[j][i - 1] + 2 * (1 - r * r) * phi[j][i] + r * r * phi[j][1] - phi[j - 1][i] - t * t * (F(j + 1, i) - F(j - 1, i)) / (phi[j + 1][i] - phi[j - 1][i]);
                    }
                } while (abs(recent - phi[j + 1][i]) > 1e-5); 
            }
            else {
                double recent = 100;
                phi[j + 1][i] = phi[j][i];
                do {
                    recent = phi[j + 1][i];
                    if (phi[j + 1][i] == phi[j - 1][i]) {
                        phi[j + 1][i] = r * r * phi[j][i - 1] + 2 * (1 - r * r) * phi[j][i] + r * r * phi[j][i + 1] - phi[j - 1][i] - t * t * (m * m * phi[j + 1][i] + g * pow(phi[j + 1][i], 3));
                    }
                    else {
                        phi[j + 1][i] = r * r * phi[j][i - 1] + 2 * (1 - r * r) * phi[j][i] + r * r * phi[j][i + 1] - phi[j - 1][i] - t * t * (F(j + 1, i) - F(j - 1, i)) / (phi[j + 1][i] - phi[j - 1][i]);
                    }
                } while (abs(recent - phi[j + 1][i]) > 1e-5); 
            }
        }
    }

    for (int j = 0; j < SIZE_T - 1; ++j) {
        for (int i = 0; i < SIZE_X-1; ++i) {
            if (i == SIZE_X - 1) {
                Energy[j] += h * (0.5 * pow((phi[j + 1][i] - phi[j][i]) / t, 2) + 0.5 * (phi[j + 1][1] - phi[j + 1][i]) / h * (phi[j][1] - phi[j][i]) / h + (F(j + 1, i) + F(j, i)) / 2);
            }
            else {
                Energy[j] += h * (0.5 * pow((phi[j + 1][i] - phi[j][i])/t, 2) + 0.5 * (phi[j + 1][i + 1] - phi[j + 1][i]) / h * (phi[j][i + 1] - phi[j][i]) / h + (F(j + 1, i) + F(j, i)) / 2);
            }
            //printf("%10.3lf, %10.3lf, %10.3lf, %10.3lf, %10.3lf, %10.3lf\n", phi[j + 1][i], phi[j][i], phi[j + 1][i + 1], phi[j][i + 1], F(j + 1, i), F(j, i));
            //printf("Energy[%d] = %lf", j, Energy[j]);
            //getchar();
        }
        Energy[SIZE_T - 1] = Energy[SIZE_T - 2];
    }
}

double eta(double x) {
    return 1. / (sqrt(2. * PI) * sigma) * exp(-(x - x_q) * (x - x_q) / (2. * sigma * sigma));
}

double theta_func(double x) {
    if (x >= 0) return 1;
    else return 0;
}

double Green_func(double x, double x_, double tm, double tm_) {
    double sum = 0;
    for (int k = -1000; k <= 1000; ++k) {
        if ((tm - tm_ - abs(x - x_ - k * L)) > 0)
            sum += _j0(m * sqrt((tm - tm_) * (tm - tm_) - (x - x_ - k * L) * (x - x_ - k * L)));
        else
            sum += 0;
    }
    //printf("res = %lf\n", sum);
    //getchar();
    return -1. / 2 * sum * al;
}

void analytical_solution() {
    for (int i = 0; i < SIZE_X; ++i) {
        for (int j = 0; j < SIZE_X; ++j) {
            for (int k = 0; k < 11; ++k) {
                phi_exact[i][k] += h * eta(x[j]) * Green_func(x[i], x[j], time_[time_moments[k]], t1);
            }
        }
        for (int k = 0; k < 11; ++k) {
            phi_asymt[i][k] = -al * time_[time_moments[k]] * eta(x[i]);
        }
    }

}

void solve_with_cond()
{
    generate_initial_conditions();
    fill_phi();
    Strauss_Vazquez();
    //Sinh_Gordon();
}