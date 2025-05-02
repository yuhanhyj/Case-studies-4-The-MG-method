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
std::vector<double> Vcycle(int N, std::vector<double> x, const std::vector<double> &b,
                           double omega, int nu, int level, int lmax);

/**
 * @brief Compute number of levels from fine grid N down to coarse_N
 */
int compute_lmax(int N, int coarse_N);

#endif // MULTIGRID_H
