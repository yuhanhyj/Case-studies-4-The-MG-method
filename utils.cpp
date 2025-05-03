#include "utils.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <cmath>

std::vector<double> matrixVectorProduct(const std::vector<double> &x, int N)
{
    int size = (N + 1) * (N + 1);
    std::vector<double> Ax(size, 0.0);
    double dx = 1.0 / N;
    double coef = 1.0 / (dx * dx);

    for (int i = 1; i < N; ++i)
    {
        for (int j = 1; j < N; ++j)
        {
            int idx = i * (N + 1) + j;
            Ax[idx] = 4.0 * x[idx] - x[idx - 1] - x[idx + 1] - x[idx - (N + 1)] - x[idx + (N + 1)];
            Ax[idx] *= coef;
        }
    }

    // Boundary condition: identity (A*x = x) at the boundary
    for (int i = 0; i <= N; ++i)
    {
        for (int j = 0; j <= N; ++j)
        {
            if (i == 0 || i == N || j == 0 || j == N)
            {
                int idx = i * (N + 1) + j;
                Ax[idx] = x[idx];
            }
        }
    }

    return Ax;
}

#include <cmath>

std::vector<double> calculateResidual(const std::vector<double> &x, const std::vector<double> &b, int N)
{
    std::vector<double> Ax = matrixVectorProduct(x, N);
    std::vector<double> r(b.size());

    for (size_t i = 0; i < b.size(); ++i)
    {
        r[i] = b[i] - Ax[i];
    }

    return r;
}

double calculateResidualNorm(const std::vector<double> &x, const std::vector<double> &b, int N)
{
    std::vector<double> r = calculateResidual(x, b, N);
    double norm = 0.0;

    for (double val : r)
    {
        norm += val * val;
    }

    return std::sqrt(norm);
}

std::vector<double> createRHSVector(int N, std::function<double(double, double)> f)
{
    double dx = 1.0 / N;
    std::vector<double> rhs((N + 1) * (N + 1), 0.0);

    for (int i = 0; i <= N; ++i)
    {
        for (int j = 0; j <= N; ++j)
        {
            int idx = i * (N + 1) + j;
            double x = j * dx;
            double y = i * dx;
            if (i != 0 && i != N && j != 0 && j != N)
            {
                rhs[idx] = f(x, y);
            }
        }
    }

    return rhs;
}

void jacobiSmoother(std::vector<double> &x, const std::vector<double> &b, int N, double omega, int iterations)
{
    int size = (N + 1) * (N + 1);
    double dx = 1.0 / N;
    double dx2 = dx * dx;
    std::vector<double> x_new(size, 0.0);

    // For small N, use a more aggressive relaxation parameter
    double adaptive_omega = omega;
    if (N <= 16)
    {
        adaptive_omega = 0.8; // More aggressive for coarse grids
    }

    for (int iter = 0; iter < iterations; ++iter)
    {
        for (int i = 1; i < N; ++i)
        {
            for (int j = 1; j < N; ++j)
            {
                int idx = i * (N + 1) + j;
                double sum = x[idx - 1] + x[idx + 1] + x[idx - (N + 1)] + x[idx + (N + 1)];
                x_new[idx] = (1.0 - adaptive_omega) * x[idx] + (adaptive_omega / 4.0) * (b[idx] * dx2 + sum);
            }
        }

        // Update x with new values
        for (int i = 1; i < N; ++i)
        {
            for (int j = 1; j < N; ++j)
            {
                int idx = i * (N + 1) + j;
                x[idx] = x_new[idx];
            }
        }
    }
}

void saveResiduals(const std::vector<double> &residuals, const std::string &filename)
{
    std::ofstream file(filename);
    for (size_t i = 0; i < residuals.size(); ++i)
    {
        file << i << " " << residuals[i] << "\n";
    }
    file.close();
    std::cout << "Saved residuals to " << filename << std::endl;
}

// ✅ 新增的 saveSolution 实现
void saveSolution(const std::vector<double> &x, int N, double dx, const std::string &filename)
{
    std::ofstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error: Cannot open file " << filename << " for writing." << std::endl;
        return;
    }

    file << N + 1 << std::endl;
    for (int i = 0; i <= N; ++i)
    {
        for (int j = 0; j <= N; ++j)
        {
            int idx = i * (N + 1) + j;
            double x_coord = j * dx;
            double y_coord = i * dx;
            file << x_coord << " " << y_coord << " " << x[idx] << "\n";
        }
    }

    file.close();
    std::cout << "Saved solution to " << filename << std::endl;
}
