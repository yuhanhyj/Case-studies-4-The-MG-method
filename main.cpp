#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <iomanip> // For table formatting
#include "multigrid.h"
#include "utils.h"

int main()
{
    // List of grid sizes to test
    std::vector<int> testNs = {16, 32, 64, 128, 256};

    // Parameters for the multigrid solver
    const double omega = 2.0 / 3.0; // Weighted Jacobi relaxation parameter
    const int nu = 3;               // Number of smoothing iterations
    const double tol = 1e-7;        // Target residual tolerance

    // Print table header for performance comparison
    std::cout << "=== Performance Comparison: 2-level vs Max-level MG (Coarsest N=8) ===" << std::endl;
    std::cout << std::left
              << std::setw(6) << "N"
              << std::setw(10) << "Method"
              << std::setw(10) << "Levels"
              << std::setw(12) << "Iterations"
              << std::setw(15) << "Coarse Solves"
              << std::setw(15) << "Runtime (ms)"
              << std::setw(15) << "Final Residual"
              << std::endl;
    std::cout << std::string(83, '-') << std::endl;

    // Loop over each grid size
    for (int N : testNs)
    {
        std::cout << "=== Solving Poisson problem for N = " << N << " ===" << std::endl;

        // Create right-hand side vector f(x) = 2π² sin(πx) sin(πy)
        std::vector<double> b = createRHSVector(N, [](double x, double y)
                                                {
            const double pi = std::acos(-1.0);
            return 2 * pi * pi * std::sin(pi * x) * std::sin(pi * y); });

        // --- 2-level Multigrid Solve ---
        std::cout << "[2-level Multigrid]" << std::endl;
        int lmax_2level = 2;
        // Initial guess x = 0
        std::vector<double> x0_2level((N + 1) * (N + 1), 0.0);
        int coarseSolves_2level = 0;
        int iterations_2level = 0;
        double runtime_2level = 0;

        std::vector<double> x_2 = multigrid_solver(
            N, x0_2level, b, omega, nu, lmax_2level, tol,
            coarseSolves_2level, iterations_2level, runtime_2level);

        double finalRes_2level = calculateResidualNorm(x_2, b, N);
        saveSolution(x_2, N, 1.0 / N, "sol_N" + std::to_string(N) + "_2level.txt");

        // --- Full Multigrid Solve (max-level) ---
        std::cout << "[Full Multigrid with coarsest level N=8]" << std::endl;
        int lmax_full = compute_lmax(N, 8); // Determine levels down to N=8
        // Fresh initial guess x = 0
        std::vector<double> x0_full((N + 1) * (N + 1), 0.0);
        int coarseSolves_full = 0;
        int iterations_full = 0;
        double runtime_full = 0;

        std::vector<double> x_full = multigrid_solver(
            N, x0_full, b, omega, nu, lmax_full, tol,
            coarseSolves_full, iterations_full, runtime_full);

        double finalRes_full = calculateResidualNorm(x_full, b, N);
        saveSolution(x_full, N, 1.0 / N, "sol_N" + std::to_string(N) + "_maxlevel.txt");

        // Print performance comparison row for 2-level
        std::cout << std::left
                  << std::setw(6) << N
                  << std::setw(10) << "2-level"
                  << std::setw(10) << lmax_2level
                  << std::setw(12) << iterations_2level
                  << std::setw(15) << coarseSolves_2level
                  << std::setw(15) << runtime_2level
                  << std::setw(15) << finalRes_2level
                  << std::endl;

        // Print performance comparison row for max-level
        std::cout << std::left
                  << std::setw(6) << N
                  << std::setw(10) << "Max-level"
                  << std::setw(10) << lmax_full
                  << std::setw(12) << iterations_full
                  << std::setw(15) << coarseSolves_full
                  << std::setw(15) << runtime_full
                  << std::setw(15) << finalRes_full
                  << std::endl;

        std::cout << std::string(83, '-') << std::endl;
    }

    std::cout << "\nDone! Run 'python3 plot_residuals.py' to visualize convergence." << std::endl;
    return 0;
}