#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

// Define the basis functions (phi_i) and their gradients
double phi(int i, double x, double y) {
    switch (i) {
        case 0: return 1.0 - x - y;
        case 1: return x;
        case 2: return y;
        default: return 0.0;
    }
}

double dphi_dx(int i, double x, double y) {
    switch (i) {
        case 0: return -1.0;
        case 1: return 1.0;
        case 2: return 0.0;
        default: return 0.0;
    }
}

double dphi_dy(int i, double x, double y) {
    switch (i) {
        case 0: return -1.0;
        case 1: return 0.0;
        case 2: return 1.0;
        default: return 0.0;
    }
}

// Define the quadrature points and weights for integration
struct QuadraturePoint {
    double x, y, weight;
};

std::vector<QuadraturePoint> quadrature_points = {
    {0.25, 0.25, 1.0},
    {0.75, 0.25, 1.0},
    {0.25, 0.75, 1.0},
    {0.75, 0.75, 1.0}
};

// Define the integrals for M, B, C1, and C2 matrices
double integrate_M(int i, int j) {
    double integral = 0.0;
    for (const auto& qp : quadrature_points) {
        integral += phi(i, qp.x, qp.y) * phi(j, qp.x, qp.y) * qp.weight;
    }
    return integral;
}

double integrate_B(int i, int j) {
    double integral = 0.0;
    for (const auto& qp : quadrature_points) {
        integral += (dphi_dx(i, qp.x, qp.y) * dphi_dx(j, qp.x, qp.y) +
                     dphi_dy(i, qp.x, qp.y) * dphi_dy(j, qp.x, qp.y)) * qp.weight;
    }
    return integral;
}

double integrate_C1(int i, int j) {
    double integral = 0.0;
    for (const auto& qp : quadrature_points) {
        integral += 0.5 * (dphi_dx(i, qp.x, qp.y) * dphi_dx(j, qp.x, qp.y) -
                           dphi_dy(i, qp.x, qp.y) * dphi_dy(j, qp.x, qp.y)) * qp.weight;
    }
    return integral;
}

double integrate_C2(int i, int j) {
    double integral = 0.0;
    for (const auto& qp : quadrature_points) {
        integral += 0.5 * (dphi_dy(i, qp.x, qp.y) * dphi_dx(j, qp.x, qp.y) -
                           dphi_dx(i, qp.x, qp.y) * dphi_dy(j, qp.x, qp.y)) * qp.weight;
    }
    return integral;
}

// Define the forcing function f(x, y)
double forcing_function(double x, double y) {
    return std::sin(M_PI * x) * std::sin(M_PI * y);
}

// Assemble the system matrix and right-hand side vector
void assemble_system(std::vector<std::vector<double>> &A, std::vector<double> &b) {
    int n = 3; // Number of basis functions

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            A[i][j] = integrate_M(i, j);
            A[i][j + n] = integrate_B(i, j);
            A[i][j + 2 * n] = 0.0;
            A[i][j + 3 * n] = 0.0;

            A[i + n][j] = integrate_B(i, j);
            A[i + n][j + n] = 0.0;
            A[i + n][j + 2 * n] = integrate_C1(i, j);
            A[i + n][j + 3 * n] = integrate_C2(i, j);

            A[i + 2 * n][j] = 0.0;
            A[i + 2 * n][j + n] = integrate_C1(i, j);
            A[i + 2 * n][j + 2 * n] = integrate_M(i, j);
            A[i + 2 * n][j + 3 * n] = 0.0;

            A[i + 3 * n][j] = 0.0;
            A[i + 3 * n][j + n] = integrate_C2(i, j);
            A[i + 3 * n][j + 2 * n] = 0.0;
            A[i + 3 * n][j + 3 * n] = integrate_M(i, j);
        }
        b[i] = 0.0;
        for (const auto& qp : quadrature_points) {
            b[i + n] += -forcing_function(qp.x, qp.y) * phi(i, qp.x, qp.y) * qp.weight;
        }
        b[i + 2 * n] = 0.0;
        b[i + 3 * n] = 0.0;
    }
}

// Solve the linear system using Gaussian elimination with partial pivoting
std::vector<double> solve_system(std::vector<std::vector<double>> &A, std::vector<double> &b) {
    int n = b.size();
    std::vector<double> x(n, 0.0);
    std::vector<int> pivot(n);

    // Initialize pivot indices
    for (int i = 0; i < n; ++i) {
        pivot[i] = i;
    }

    // Gaussian elimination with partial pivoting
    for (int i = 0; i < n; ++i) {
        // Find pivot
        int max_row = i;
        for (int k = i + 1; k < n; ++k) {
            if (std::abs(A[pivot[k]][i]) > std::abs(A[pivot[max_row]][i])) {
                max_row = k;
            }
        }
        std::swap(pivot[i], pivot[max_row]);

        if (A[pivot[i]][i] == 0.0) {
            std::cerr << "Zero pivot encountered at row " << i << std::endl;
            exit(1);
        }

        for (int j = i + 1; j < n; ++j) {
            double factor = A[pivot[j]][i] / A[pivot[i]][i];
            for (int k = i; k < n; ++k) {
                A[pivot[j]][k] -= factor * A[pivot[i]][k];
            }
            b[pivot[j]] -= factor * b[pivot[i]];
        }
    }

    // Back substitution
    for (int i = n - 1; i >= 0; --i) {
        x[i] = b[pivot[i]] / A[pivot[i]][i];
        for (int j = i - 1; j >= 0; --j) {
            b[pivot[j]] -= A[pivot[j]][i] * x[i];
        }
    }

    return x;
}

int main() {
    int n = 3; // Number of basis functions
    std::vector<std::vector<double>> A(4 * n, std::vector<double>(4 * n, 0.0));
    std::vector<double> b(4 * n, 0.0);

    // Assemble the system
    assemble_system(A, b);

    // Solve the system
    std::vector<double> solution = solve_system(A, b);

    // Output the solution
    std::cout << "Solution:" << std::endl;
    for (double val : solution) {
        std::cout << val << " ";
    }
    std::cout << std::endl;

    return 0;
}
