#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <functional>

/**
 * @brief Matrix-vector product using 5-point stencil for 2D Poisson problem.
 *        No explicit matrix is stored.
 * @param x Input vector representing grid values
 * @param N Grid size (number of intervals in each direction)
 * @return A * x result as a vector
 */
std::vector<double> matrixVectorProduct(const std::vector<double> &x, int N);

/**
 * @brief Generate the right-hand side vector f(x, y) on a 2D grid.
 * @param N Grid size (number of intervals in each direction)
 * @param f Function f(x, y) for RHS of Poisson equation
 * @return Discretized RHS vector
 */
std::vector<double> createRHSVector(int N, std::function<double(double, double)> f);

/**
 * @brief Perform weighted Jacobi smoothing.
 * @param x Solution vector to be updated
 * @param b Right-hand side vector
 * @param N Grid size
 * @param omega Relaxation parameter
 * @param iterations Number of smoothing iterations
 */
void jacobiSmoother(std::vector<double> &x, const std::vector<double> &b, int N, double omega, int iterations);

#endif // UTILS_H
