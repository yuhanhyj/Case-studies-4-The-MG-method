#ifndef MULTIGRID_H
#define MULTIGRID_H

#include <vector>

/**
 * @brief Solve Poisson problem using multigrid with V-cycle.
 * @param N Grid size
 * @param x Initial guess
 * @param b Right-hand side vector
 * @param omega Weighted Jacobi relaxation parameter
 * @param nu Number of smoothing iterations
 * @param lmax Number of grid levels
 * @param tol Target residual tolerance
 * @return Solution vector
 */
std::vector<double> multigrid_solver(int N, std::vector<double> x, const std::vector<double> &b, double omega, int nu, int lmax, double tol);

/**
 * @brief One recursive V-cycle
 */
std::vector<double> Vcycle(int N, std::vector<double> x, const std::vector<double> &b, double omega, int nu, int level, int lmax);

/**
 * @brief Compute number of levels from fine grid N down to coarse_N
 */
int compute_lmax(int N, int coarse_N);

/**
 * @brief Restrict residual from fine grid (N) to coarse grid (N/2)
 * @param r_fine Residual vector on fine grid
 * @param N_fine Fine grid size
 * @return Restricted residual on coarse grid
 */
std::vector<double> restrictResidual(const std::vector<double> &r_fine, int N_fine);

/**
 * @brief Prolong correction from coarse grid (N) to fine grid (2N)
 * @param e_coarse Error correction on coarse grid
 * @param N_coarse Coarse grid size
 * @return Prolongated correction on fine grid
 */
std::vector<double> prolongCorrection(const std::vector<double> &e_coarse, int N_coarse);

#endif // MULTIGRID_H
