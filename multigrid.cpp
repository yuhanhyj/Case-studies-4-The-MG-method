#include "multigrid.h"
#include "utils.h"
#include <cmath>
#include <iostream>

// Solve using repeated V-cycles
std::vector<double> multigrid_solver(int N, std::vector<double> x, const std::vector<double> &b, double omega, int nu, int lmax, double tol)
{
    const int maxIters = 50;

    for (int iter = 0; iter < maxIters; ++iter)
    {
        x = Vcycle(N, x, b, omega, nu, 1, lmax);

        double resNorm = calculateResidualNorm(x, b, N);
        std::cout << "  Iteration " << iter + 1 << ", residual = " << resNorm << std::endl;

        if (resNorm < tol)
            break;
    }

    return x;
}

// Initial V-cycle structure (we'll add restriction/prolong later)
std::vector<double> Vcycle(int N, std::vector<double> x, const std::vector<double> &b,、double omega, int nu, int level, int lmax)
{
    // Step 1: smoothing (not yet implemented)
    // smooth(N, x, b, omega, nu);
    // Pre-smoothing
    jacobiSmoother(x, b, N, omega, nu);

    // Step 2–5: placeholder, returns x directly for now
    return x;
}

// Compute how many levels from N down to coarse_N
int compute_lmax(int N, int coarse_N)
{
    int levels = 1;
    while ((N >>= 1) >= coarse_N)
        ++levels;
    return levels;
}
