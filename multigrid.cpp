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

// Restriction: Full-weighting from fine (N) to coarse (N/2)
std::vector<double> restrictResidual(const std::vector<double> &r_fine, int N_fine)
{
    int N_coarse = N_fine / 2;
    std::vector<double> r_coarse((N_coarse + 1) * (N_coarse + 1), 0.0);

    for (int i = 1; i < N_coarse; ++i)
    {
        for (int j = 1; j < N_coarse; ++j)
        {
            int ic = i * (N_coarse + 1) + j;

            int fi = 2 * i;
            int fj = 2 * j;
            int fc = fi * (N_fine + 1) + fj;

            r_coarse[ic] = 0.25 * r_fine[fc] + 0.125 * (r_fine[fc - 1] + r_fine[fc + 1] + r_fine[fc - (N_fine + 1)] + r_fine[fc + (N_fine + 1)]) + 0.0625 * (r_fine[fc - (N_fine + 2)] + r_fine[fc - (N_fine)] + r_fine[fc + (N_fine + 2)] + r_fine[fc + (N_fine)]);
        }
    }

    return r_coarse;
}

// Prolongation: Bilinear interpolation from coarse (N) to fine (2N)
std::vector<double> prolongCorrection(const std::vector<double> &e_coarse, int N_coarse)
{
    int N_fine = 2 * N_coarse;
    std::vector<double> e_fine((N_fine + 1) * (N_fine + 1), 0.0);

    for (int i = 0; i <= N_coarse; ++i)
    {
        for (int j = 0; j <= N_coarse; ++j)
        {
            int coarse_idx = i * (N_coarse + 1) + j;
            int fi = 2 * i;
            int fj = 2 * j;
            int fine_idx = fi * (N_fine + 1) + fj;
            e_fine[fine_idx] = e_coarse[coarse_idx];
        }
    }

    // Interpolate horizontally and vertically
    for (int i = 1; i < N_fine; i += 2)
    {
        for (int j = 0; j <= N_fine; ++j)
        {
            int idx = i * (N_fine + 1) + j;
            e_fine[idx] = 0.5 * (e_fine[(i - 1) * (N_fine + 1) + j] + e_fine[(i + 1) * (N_fine + 1) + j]);
        }
    }
    for (int i = 0; i <= N_fine; ++i)
    {
        for (int j = 1; j < N_fine; j += 2)
        {
            int idx = i * (N_fine + 1) + j;
            e_fine[idx] = 0.5 * (e_fine[i * (N_fine + 1) + (j - 1)] + e_fine[i * (N_fine + 1) + (j + 1)]);
        }
    }

    return e_fine;
}
