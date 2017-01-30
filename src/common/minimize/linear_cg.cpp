#include "common/minimize/minimize.hpp"
#include "common/util/blas.hpp"

#include <iostream>
#include <cmath>

using namespace std;
using namespace sp2;

std::vector<double> minimize::linear_cg(vector_fn_t matrix_fn,
     std::vector<double> b, double tolerance)
{
    if (tolerance <= 0)
        throw invalid_argument("tolerance must be positive");

    // dimension
    auto n = b.size();

    // current solution x, residual, and basis vector
    auto x = vector<double>(n, 0), // x_0 = 0
        residual = b,              // r_0 = b - Ax_0 = b
        basis = residual;          // p_0 = r_0

    // residual dot product with itself
    auto residual_dot = vdot(residual, residual);

    while (max_norm(residual) > tolerance)
    {
        // Ap_k
        auto basis_A = matrix_fn(basis);

        // alpha_k = (r_k^T . r_k) / (p_k^T . Ap_k)
        auto alpha = residual_dot / vdot(basis, basis_A);

        ////////
        // get the new x and residual

        // x_{k + 1} = x_k + alpha_k * p_k
        vaxpy(alpha, basis, x);

        // r_{k + 1} = r_k + alpha_k * Ap_k
        vaxpy(-alpha, basis_A, residual);

        ////////
        // get beta and calculate the new basis

        // beta_k = (r_{k + 1}^T . r_{k + 1}) / (r_k^T . r_k)
        auto new_res_dot = vdot(residual, residual),
            beta = new_res_dot / residual_dot;

        residual_dot = new_res_dot;

        // p_{k + 1} = r_{k + 1} + beta_k * p_k
        auto new_basis = residual;
        vaxpy(beta, basis, new_basis);

        basis = new_basis;
    }

    return x;
}
