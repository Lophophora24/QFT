#include "MonteCarlo.h"

//----------------- Что будем усреднять --------------------//
// Корреляционные функции (L/2, x)
double Corr_function_t_0[SIZE_X];
double Corr_function_t_1[SIZE_X];
double Corr_function_t_2[SIZE_X];
double Corr_function_t_3[SIZE_X];
double Corr_function_t_4[SIZE_X];
double Corr_function_t_5[SIZE_X];

// Корреляционные функции (x, x)
double Corr_function_t_0_0[SIZE_X];
double Corr_function_t_1_1[SIZE_X];
double Corr_function_t_2_2[SIZE_X];
double Corr_function_t_3_3[SIZE_X];
double Corr_function_t_4_4[SIZE_X];
double Corr_function_t_5_5[SIZE_X];

// Корреляционные функции (0, x)
double Corr_function_t_0_0_0[SIZE_X];
double Corr_function_t_1_1_1[SIZE_X];
double Corr_function_t_2_2_2[SIZE_X];
double Corr_function_t_3_3_3[SIZE_X];
double Corr_function_t_4_4_4[SIZE_X];
double Corr_function_t_5_5_5[SIZE_X];
double Result = 0;

double phi_aver[SIZE_X][11];
double energy_aver[SIZE_T];

double phi_4_aver[SIZE_T];

void average()
{
	//FILE* samples;
	//fopen_s(&samples, "samples_new.txt", "a+");
	printf("Start to solve equations\n");
	for (int i = 0; i < M; ++i) {
		printf("...\n");
		solve_with_cond();  // получили phi с гауссовыми начальными условиями
		printf("Solved equation number # %d\n", i);
		
		/*fprintf(samples, "Sample #%d\n", i);
		for (int i = 0; i < SIZE_T; ++i) {
			for (int j = 0; j < SIZE_X; ++j) {
				fprintf(samples, "%10.3lf", phi[i][j]);
			}
			fprintf(samples, "\n");
		}*/
		

		// Wigner_func();   // посчитать Вигнеровский функционал на начальных условиях
		                    // (функции phi_0 и pi_0)

		// усреднить единицу и убедиться, что усреднение дает единицу

		//Result += Wigner_func() / M;

		//Result += (phi[0][0] * phi[0][0]) / M;
		
		//здесь считаются корреляционные функции (L/2,x), (x, x), (0, x)
		/*
		for (int j = 0; j < SIZE_X; ++j) {
			Corr_function_t_0[j] += (phi[0][SIZE_X / 2] * phi[0][j]) / M;
			Corr_function_t_1[j] += (phi[200][SIZE_X / 2] * phi[200][j]) / M;
			Corr_function_t_2[j] += (phi[400][SIZE_X / 2] * phi[400][j]) / M;
			Corr_function_t_3[j] += (phi[600][SIZE_X / 2] * phi[600][j]) / M;
			Corr_function_t_4[j] += (phi[800][SIZE_X / 2] * phi[800][j]) / M;
			Corr_function_t_5[j] += (phi[900][SIZE_X / 2] * phi[900][j]) / M;

			Corr_function_t_0_0[j] += (phi[0][j] * phi[0][j]) / M;
			Corr_function_t_1_1[j] += (phi[200][j] * phi[200][j]) / M;
			Corr_function_t_2_2[j] += (phi[400][j] * phi[400][j]) / M;
			Corr_function_t_3_3[j] += (phi[600][j] * phi[600][j]) / M;
			Corr_function_t_4_4[j] += (phi[800][j] * phi[800][j]) / M;
			Corr_function_t_5_5[j] += (phi[900][j] * phi[900][j]) / M;

			Corr_function_t_0_0_0[j] += (phi[0][0] * phi[0][j]) / M;
			Corr_function_t_1_1_1[j] += (phi[200][0] * phi[200][j]) / M;
			Corr_function_t_2_2_2[j] += (phi[400][0] * phi[400][j]) / M;
			Corr_function_t_3_3_3[j] += (phi[600][0] * phi[600][j]) / M;
			Corr_function_t_4_4_4[j] += (phi[800][0] * phi[800][j]) / M;
			Corr_function_t_5_5_5[j] += (phi[900][0] * phi[900][j]) / M;
		}
		*/

		for (int j = 0; j < SIZE_T; ++j) {
			energy_aver[j] += Energy[j] / M;

			//----------------------------------------------------------

			double int_multiply = 0;
			for (int i = 0; i < SIZE_X; ++i) {
				int_multiply += g_0 * pow(phi[j][i], 4) * h / (4. * TMP);
			}

			energy_aver[j] *= exp(-int_multiply);

			//----------------------------------------------------------
		}

		//------------------------------------------------------
		
		for (int j = 0; j < SIZE_T; ++j) {
			double int_multiply = 0;
			for (int i = 0; i < SIZE_X; ++i) {
				int_multiply += g_0 * pow(phi[j][i], 4) * h / (4. * TMP);
			}

			phi_4_aver[j] += exp(-int_multiply) / M;
		}

		for (int j = 0; j < SIZE_T; ++j) {
			energy_aver[j] /= phi_4_aver[j];
		}

		//------------------------------------------------------
		
		for (int k = 0; k < 11; ++k) {
			//--------------------------------------------------------------------

			double multip = 0;

			for (int j = 0; j < SIZE_X; ++j) {
				multip += g_0 * pow(phi[time_moments[k]][j], 4) * h / (4. * TMP);
			}



			//--------------------------------------------------------------------
			for (int j = 0; j < SIZE_X - 1; ++j) {
				//phi_aver[j][k] += (phi[time_moments[k]][j] * exp(phi_averaged_with_eta)) / M;
				/*
				phi_aver[j][k] += (phi[time_moments[k]][j]) / M;
				phi_aver[SIZE_X - 1][k] = phi_aver[0][k];
				*/
				//if (j == SIZE_X - 1) {
					//phi_aver[j][k] += (0.5 * pow((phi[time_moments[k] + 1][j] - phi[time_moments[k]][j]) / t, 2) + 0.5 * (phi[time_moments[k] + 1][1] - phi[time_moments[k] + 1][j]) / h * (phi[time_moments[k]][1] - phi[time_moments[k]][j]) / h + (F(time_moments[k] + 1, j) + F(time_moments[k], j)) / 2) / M;
					//phi_aver[j][k] += 0.5 * (0.5 * pow((phi[time_moments[k] + 1][j] - phi[time_moments[k]][j]) / t, 2) + 0.5 * (phi[time_moments[k] + 1][1] - phi[time_moments[k] + 1][j]) / h * (phi[time_moments[k]][1] - phi[time_moments[k]][j]) / h + (F(time_moments[k] + 1, j) + F(time_moments[k], j)) / 2) / M + 0.5 * (0.5 * pow((phi[time_moments[k] + 1][j] - phi[time_moments[k]][j]) / t, 2) + 0.5 * (phi[time_moments[k] + 1][j - 1] - phi[time_moments[k] + 1][j]) / h * (phi[time_moments[k]][j - 1] - phi[time_moments[k]][j]) / h + (F(time_moments[k] + 1, j) + F(time_moments[k], j)) / 2) / M;
				//}
				
				if (j == 0) {
					phi_aver[j][k] += 0.5 * (0.5 * pow((phi[time_moments[k] + 1][j] - phi[time_moments[k]][j]) / t, 2) + 0.5 * (phi[time_moments[k] + 1][j + 1] - phi[time_moments[k] + 1][j]) / h * (phi[time_moments[k]][j + 1] - phi[time_moments[k]][j]) / h + (F(time_moments[k] + 1, j) + F(time_moments[k], j)) / 2) / M + 0.5 * (0.5 * pow((phi[time_moments[k] + 1][j] - phi[time_moments[k]][j]) / t, 2) + 0.5 * (phi[time_moments[k] + 1][SIZE_X-2] - phi[time_moments[k] + 1][j]) / h * (phi[time_moments[k]][SIZE_X-2] - phi[time_moments[k]][j]) / h + (F(time_moments[k] + 1, j) + F(time_moments[k], j)) / 2) * exp(-multip) / M;
				}
				else {
					phi_aver[j][k] += 0.5 * (0.5 * pow((phi[time_moments[k] + 1][j] - phi[time_moments[k]][j]) / t, 2) + 0.5 * (phi[time_moments[k] + 1][j + 1] - phi[time_moments[k] + 1][j]) / h * (phi[time_moments[k]][j + 1] - phi[time_moments[k]][j]) / h + (F(time_moments[k] + 1, j) + F(time_moments[k], j)) / 2) / M + 0.5 * (0.5 * pow((phi[time_moments[k] + 1][j] - phi[time_moments[k]][j]) / t, 2) + 0.5 * (phi[time_moments[k] + 1][j - 1] - phi[time_moments[k] + 1][j]) / h * (phi[time_moments[k]][j - 1] - phi[time_moments[k]][j]) / h + (F(time_moments[k] + 1, j) + F(time_moments[k], j)) / 2) * exp(-multip) / M;	
				}

				phi_aver[SIZE_X - 1][k] = phi_aver[0][k];
				

				/*if (j == SIZE_X - 1) {
					phi_aver[j][k] += (0.5 * pow((phi[time_moments[k] + 1][j] - phi[time_moments[k]][j]) / t, 2) + 0.5 * pow((phi[time_moments[k]][1] - phi[time_moments[k]][j]) / h, 2) + 0.5 * F(time_moments[k], j)) / M;
				}
				else {
					phi_aver[j][k] += (0.5 * pow((phi[time_moments[k] + 1][j] - phi[time_moments[k]][j]) / t, 2) + 0.5 * pow((phi[time_moments[k]][j + 1] - phi[time_moments[k]][j]) / h, 2) + 0.5 * F(time_moments[k], j)) / M;
				}*/

				//phi_aver[j][k] *= exp(-multip);
			}
			//phi_aver[SIZE_X-1][k] *= exp(-multip);
		}
		
		// сделать усреднение <phi_n, phi_m> это должно совпадать с коррел. матрицей???
		// сделать усреднение <phi_i, phi_j, phi_n, phi_m> = <phi_i, phi_j><phi_n, phi_m> +
		// + <><> + <><> проверить, что это соотношение выполняется 
		//std::cout << "i = " << i << "  res = " << Result << '\n';
	}

	for (int k = 0; k < 11; ++k) {
		for (int j = 0; j < SIZE_X - 1; ++j) {
			phi_aver[j][k] /= phi_4_aver[time_moments[k]];
		}
	}

	//fclose(samples);

	//std::cout << "Result: " << Result << '\n';
}