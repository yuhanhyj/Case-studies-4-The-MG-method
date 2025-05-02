#include "utils.h"
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
