#include "NdimGauss.h"

double phi_0[SIZE_X];
double pi_0[SIZE_X];

double phi_averaged_with_eta = 0;

Eigen::VectorXd mean_phi(SIZE_X-1);             // вектор средних для поля phi
Eigen::MatrixXd covar_phi(SIZE_X-1, SIZE_X-1);	  // матрица ковариаций для поля phi
Eigen::VectorXd mean_pi(SIZE_X-1);              // вектор средних для поля pi
Eigen::MatrixXd covar_pi(SIZE_X-1, SIZE_X-1);	  // матрица ковариаций для поля pi

Eigen::MatrixXd A(SIZE_X-1, SIZE_X-1);            // матрицы для вигнеровского функционала
Eigen::MatrixXd B(SIZE_X-1, SIZE_X-1);   

void init_()
{
	for (int i = 0; i < SIZE_X-1; ++i) {
		mean_phi(i) = 0;
		mean_pi(i) = 0;
		for (int j = 0; j < SIZE_X-1; ++j) {
			covar_phi(i, j) = 0;
			covar_pi(i, j) = 0;
		}
	}

	for (int i = 0; i < SIZE_X-1; ++i) {
		for (int j = 0; j < SIZE_X-1; ++j) {
			A(i, j) = 0;
			B(i, j) = 0;
			if (i == j) {
				A(i, j) = h / TMP;                        // здесь заменил - на +, хотя должен быть - (в начале формулы), а может и не должен
				B(i, j) = (m * m * h) / TMP + 2. / (h * TMP);   // и здесь

				//A(i, j) = 1;
				//B(i, j) = 1;

				//if ((i == 0) || (i == (SIZE_X - 1))) {
					//B(i, j) += (-1. / (h * TMP));
				//}
			}
			if (abs(i - j) == 1) {
				B(i, j) = -1. / (TMP * h);         // и здесь тоже
			}
		}
	}

	B(0, SIZE_X-2) = (-1. / (h * TMP));
	B(SIZE_X-2, 0) = (-1. / (h * TMP));

	covar_pi = A.inverse();
	covar_phi = B.inverse();
/*
	//std::cout << covar_phi(0, 0);
	FILE* phi_covar;
	fopen_s(&phi_covar, "covar_phi.txt", "w+");
	
	if (phi_covar != 0) {
		for (int i = 0; i < SIZE_X; ++i) {
			for (int j = 0; j < SIZE_X; ++j) {
				fprintf(phi_covar, "%10.3lf", covar_phi(i, j));
			}
			fprintf(phi_covar, "\n");
		}
		fclose(phi_covar);
	}
	else {
		printf("phi_covar = 0");
	}
*/
}

int nn = 1;     // How many samples (columns) to draw

Eigen::MatrixXd normTransform_phi(SIZE_X-1, SIZE_X-1);
Eigen::MatrixXd normTransform_pi(SIZE_X-1, SIZE_X-1);
Eigen::LLT<Eigen::MatrixXd> cholSolver_phi(covar_phi);
Eigen::LLT<Eigen::MatrixXd> cholSolver_pi(covar_pi);

void init_for_generating_initial_conditions() {

	
	//------------------------Разложение для поля phi------------------------//
	if (cholSolver_phi.info() == Eigen::Success) {
		// Use cholesky solver
		normTransform_phi = cholSolver_phi.matrixL();
	}
	else {
		// Use eigen solver
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver_phi(covar_phi);
		normTransform_phi = eigenSolver_phi.eigenvectors()
			* eigenSolver_phi.eigenvalues().cwiseSqrt().asDiagonal();
	}

	//------------------------Разложение для поля pi-------------------------//
	if (cholSolver_pi.info() == Eigen::Success) {
		// Use cholesky solver
		normTransform_pi = cholSolver_pi.matrixL();
	}
	else {
		// Use eigen solver
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver_pi(covar_pi);
		normTransform_pi = eigenSolver_pi.eigenvectors()
			* eigenSolver_pi.eigenvalues().cwiseSqrt().asDiagonal();
	}
}

void generate_initial_conditions() {
	// phi_0 = ...
	// pi_0 = ...
	phi_averaged_with_eta = 0;

	//for (int i = 0; i < SIZE_X; ++i) {
	//	phi_0[i] = (1 + cos(2 * PI * x[i] / (x2 - x1))) * 0.1;
	//	pi_0[i] = 0;
	//}


	// We can only use the cholesky decomposition if
	// the covariance matrix is symmetric, pos-definite.
	// But a covariance matrix might be pos-semi-definite.
	// In that case, we'll go to an EigenSolver


	Eigen::internal::scalar_normal_dist_op<double> randN; // Gaussian functor
	Eigen::internal::scalar_normal_dist_op<double>::rng.seed(std::chrono::system_clock::now().time_since_epoch().count()); // Seed the rng

	Eigen::MatrixXd samples_phi = (normTransform_phi
		* Eigen::MatrixXd::NullaryExpr(SIZE_X-1, nn, randN)).colwise()
		+ mean_phi;

	Eigen::MatrixXd samples_pi = (normTransform_pi
		* Eigen::MatrixXd::NullaryExpr(SIZE_X-1, nn, randN)).colwise()
		+ mean_pi;

	//std::cout << "Samples_phi\n" << samples_phi << std::endl;

	//std::cout << "Samples_pi\n" << samples_pi << std::endl;
	
	for (int i = 0; i < SIZE_X-1; ++i) {
		phi_averaged_with_eta += eta(x[i]) * samples_phi(i) * h;
	}
	
	for (int i = 0; i < SIZE_X-1; ++i) {
		phi_0[i] = samples_phi(i);
		//phi_0[i] = 0.1 * (1 + cos(2 * PI / 1.28 * x[i]));
		//printf("Phi_averaged = %f\n", phi_averaged_with_eta);
		pi_0[i] = samples_pi(i) - 1. * al * eta(x[i]) * pow(phi_averaged_with_eta, 0);
		//pi_0[i] = samples_pi(i);
		//pi_0[i] = 0;
		//pi_0[i] = samples_pi(i) - 1. * al * eta(x[i]);

		//std::cout << phi_0[i] << '\n';
	}

	//phi_averaged_with_eta = 0;

	phi_0[SIZE_X - 1] = phi_0[0];
	pi_0[SIZE_X - 1] = pi_0[0];
}







