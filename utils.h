#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <string>
#include <functional>

/**
 * @brief Matrix-vector product using 5-point stencil for 2D Poisson problem.
 */
std::vector<double> matrixVectorProduct(const std::vector<double> &x, int N);

/**
 * @brief Compute residual vector r = b - A * x
 */
std::vector<double> calculateResidual(const std::vector<double> &x, const std::vector<double> &b, int N);

/**
 * @brief Compute L2 norm of residual vector
 */
double calculateResidualNorm(const std::vector<double> &x, const std::vector<double> &b, int N);

/**
 * @brief Generate the right-hand side vector f(x, y) on a 2D grid.
 */
std::vector<double> createRHSVector(int N, std::function<double(double, double)> f);

/**
 * @brief Perform weighted Jacobi smoothing.
 */
void jacobiSmoother(std::vector<double> &x, const std::vector<double> &b, int N, double omega, int iterations);

/**
 * @brief Save solution vector to file in (x, y, u) format.
 */
void saveSolution(const std::vector<double> &x, int N, double dx, const std::string &filename);

/**
 * @brief Save residuals to file for convergence plot
 */
void saveResiduals(const std::vector<double> &residuals, const std::string &filename);

#endif // UTILS_H
