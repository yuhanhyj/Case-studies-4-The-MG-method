#include "multigrid.h"
#include "utils.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <chrono>

// Solve using repeated V-cycles
std::vector<double> multigrid_solver(int N, std::vector<double> x, const std::vector<double> &b, double omega, int nu, int lmax, double tol, int &coarseSolves, int &iterations, double &runtime)
{
    const int maxIters = 100; // Increased from 50 to give more chances to converge
    std::vector<double> residuals;
    coarseSolves = 0;
    iterations = 0;

    // Track stagnation and divergence
    double prevResNorm = std::numeric_limits<double>::max();
    int stagnationCount = 0;
    const double stagnationTolerance = 0.05; // 5% improvement threshold
    const int maxStagnationIters = 5;        // Max consecutive stagnation iterations

    // Start timer
    auto startTime = std::chrono::high_resolution_clock::now();

    for (int iter = 0; iter < maxIters; ++iter)
    {
        x = Vcycle(N, x, b, omega, nu, 1, lmax, coarseSolves);

        double resNorm = calculateResidualNorm(x, b, N);
        residuals.push_back(resNorm);
        iterations = iter + 1;

        std::cout << "  Iteration " << iter + 1 << ", residual = " << resNorm << std::endl;

        // Check convergence
        if (resNorm < tol)
        {
            std::cout << "  Converged to tolerance " << tol << std::endl;
            break;
        }

        // Check for divergence
        if (iter > 0 && resNorm > prevResNorm * 1.1)
        { // 10% increase in residual
            std::cout << "  Divergence detected! Stopping iterations." << std::endl;
            break;
        }

        // Check for stagnation
        if (iter > 0 && prevResNorm - resNorm < prevResNorm * stagnationTolerance)
        {
            stagnationCount++;
            if (stagnationCount >= maxStagnationIters)
            {
                std::cout << "  Stagnation detected! Insufficient progress over "
                          << maxStagnationIters << " iterations." << std::endl;
                break;
            }
        }
        else
        {
            stagnationCount = 0; // Reset counter if good progress
        }

        prevResNorm = resNorm;
    }

    // Calculate runtime
    auto endTime = std::chrono::high_resolution_clock::now();
    runtime = std::chrono::duration<double, std::milli>(endTime - startTime).count();

    // Report summary
    std::cout << "  Total iterations: " << iterations << std::endl;
    std::cout << "  Total coarse grid solves: " << coarseSolves << std::endl;
    std::cout << "  Runtime: " << runtime << " ms" << std::endl;

    // Save residuals to file
    std::string filename = "residual_N" + std::to_string(N) + "_lmax" + std::to_string(lmax);
    // Add a method identifier to distinguish between 2-level and max-level
    if (lmax == 2)
    {
        filename += "_2level.txt";
    }
    else
    {
        filename += "_maxlevel.txt";
    }
    saveResiduals(residuals, filename);

    return x;
}

// Initial V-cycle structure with coarse solve counting
std::vector<double> Vcycle(int N, std::vector<double> x, const std::vector<double> &b, double omega, int nu, int level, int lmax, int &coarseSolves)
{
    // 1. Pre-smoothing
    jacobiSmoother(x, b, N, omega, nu);

    // 2. Compute residual r = b - A * x
    std::vector<double> r = calculateResidual(x, b, N);

    // 3. Restrict residual to coarse grid
    std::vector<double> r_coarse = restrictResidual(r, N);
    int N_coarse = N / 2;

    // 4. Coarse grid solve (recursively or directly)
    std::vector<double> e_coarse;
    if (level + 1 == lmax)
    {
        // Coarsest level: direct solve using Jacobi
        std::vector<double> zero((N_coarse + 1) * (N_coarse + 1), 0.0);
        e_coarse = zero;
        // Use more iterations for better coarse grid solve
        jacobiSmoother(e_coarse, r_coarse, N_coarse, omega, 100);
        coarseSolves++; // Count this coarse solve
    }
    else
    {
        std::vector<double> zero((N_coarse + 1) * (N_coarse + 1), 0.0);
        e_coarse = Vcycle(N_coarse, zero, r_coarse, omega, nu, level + 1, lmax, coarseSolves);
    }

    // 5. Prolongate error and correct fine grid solution
    std::vector<double> e_fine = prolongCorrection(e_coarse, N_coarse);
    for (size_t i = 0; i < x.size(); ++i)
    {
        x[i] += e_fine[i];
    }

    // 6. Post-smoothing
    jacobiSmoother(x, b, N, omega, nu);

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