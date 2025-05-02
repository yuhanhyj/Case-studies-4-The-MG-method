// src/main.cpp

#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
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

    // Loop over each grid size
    for (int N : testNs)
    {
        std::cout << "=== Solving Poisson problem for N = " << N << " ===" << std::endl;

        // Create right-hand side vector f(x) = 2π² sin(πx) sin(πy)
        std::vector<double> b = createRHSVector(N, [](double x, double y)
                                                {
            const double pi = std::acos(-1.0);
            return 2 * pi * pi * std::sin(pi * x) * std::sin(pi * y); });

        // Initial guess x = 0
        std::vector<double> x0((N + 1) * (N + 1), 0.0);

        // --- 2-level Multigrid Solve ---
        std::cout << "[2-level Multigrid]" << std::endl;
        int lmax_2level = 2;
        std::vector<double> x_2 = multigrid_solver(N, x0, b, omega, nu, lmax_2level, tol);
        saveSolution(x_2, N, 1.0 / N, "results/sol_N" + std::to_string(N) + "_2level.txt");

        // --- Full Multigrid Solve (max-level) ---
        std::cout << "[Full Multigrid with coarsest level N=8]" << std::endl;
        int lmax_full = compute_lmax(N, 8); // Determine levels down to N=8
        std::vector<double> x_full = multigrid_solver(N, x0, b, omega, nu, lmax_full, tol);
        saveSolution(x_full, N, 1.0 / N, "results/sol_N" + std::to_string(N) + "_maxlevel.txt");

        std::cout << std::endl;
    }

    return 0;
}
