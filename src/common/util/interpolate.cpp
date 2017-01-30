#include <iostream>
#include <utility>
#include <vector>

#include <stdexcept>

#include "common/util/math.hpp"
#include "common/util/templates.hpp"

// for k = 1 ... min(m,n):
//   Find the k-th pivot:
//   i_max  := argmax (i = k ... m, abs(A[i, k]))
//   if A[i_max, k] = 0
//     error "Matrix is singular!"
//   swap rows(k, i_max)
//   Do for all rows below pivot:
//   for i = k + 1 ... m:
//     f := A[i, k] / A[k, k]
//     Do for all remaining elements in current row:
//     for j = k + 1 ... n:
//       A[i, j]  := A[i, j] - A[k, j] * f
//     Fill lower triangular matrix with zeros:
//     A[i, k]  := 0

// TODO: Finish gauss_elim
void gauss_elim(const std::size_t m, const std::size_t n,
    std::vector<double> &matrix)
{
    const auto min_dim = std::min(m, n);
    if (min_dim == 0)
        return;

    // turn matrix into row echelon form
    for (std::size_t i = 0; i < std::min(m, n); ++i)
    {
        // find pivot
        auto max_elem = std::make_pair(matrix[i], std::size_t(0));
        for (auto j = i; j < m; ++j)
            max_elem = std::max(max_elem,
                std::make_pair(matrix[j * n + i], j));

        if (max_elem.first == 0)
            throw std::runtime_error("Matrix is singular.");

        // swap rows
        using std::swap;
        if (i != max_elem.first)
            for (std::size_t j = 0; j < n; ++j)
                swap(matrix[max_elem.second * n + j], matrix[i * n + j]);

        for (auto j = i + 1; j < m; ++j)
        {
            auto multiplier = matrix[j * n + i] / matrix[i * n + i];
            for (auto k = i + 1; k < n; ++k)
                matrix[j * n + k] -= matrix[i * n + k] * multiplier;

            matrix[j * n + i] = 0;
        }
    }

//    // do back-substitution to turn into reduced row echelon
//    for (std::size_t i = m - 1; ; --i)
//    {
//        for (std::size_t j = 0;
//
//        if (i == 0)
//            break;
//    }
}

struct interp_point_t
{
    /// the data value at the specified coordinate
    double value;

    /// data vector vector of length (dimension):
    ///     ((x_1, d_1), (x_2, d_2), ... , (x_n, d_n))
    /// consists of pairs of (T, int) indicating the coordinate value
    /// and the derivative order (d_n) for the coordinate
    std::vector<std::pair<double, unsigned int>> data;
};

/// input data class for get_interp_coeff()
struct interp_data_t
{
    std::size_t dimension,
        degree;

    std::vector<interp_point_t> data;

    std::vector<double> get_coefficients() const
    {
        // number of terms in the polynomial
        const auto n_terms = sp2::int_pow(degree + 1, dimension);

        // matrix that will be used to calculate the polynomial
        // coefficients for the given interpolation data
        std::vector<double> interp_matrix((n_terms + 1) * (n_terms + 1), 0);

        // for each data point, we generate a row of the interpolation matrix
        for (std::size_t i = 0; i < data.size(); ++i)
        {
            // first set the value for the data point (last value in the row)
            interp_matrix[i * (n_terms + 1) + n_terms] = data[i].value;

            // matrix row where data will be output
            double *output = &interp_matrix[i * (n_terms + 1)];

            std::size_t output_len = 0;
            // done in reverse due to something to do with outer products...
            // can't remember why
            for (auto &pair : sp2::reverse_range(data[i].data))
            {
                const auto coord = pair.first;
                const auto deriv = pair.second;

                // essentially get the 'coord' components of the polynomial
                // terms, factorial coefficients are due to derivatives
                std::vector<double> comp(degree + 1, 0);

                double power = 1;
                for (auto k = deriv; k < comp.size(); ++k, power *= coord)
                    comp[k] = sp2::factorial_diff(k, k - deriv) * power;

                if (output_len == 0)
                {
                    // first iteration
                    std::copy(comp.begin(), comp.end(), output);
                    output_len = comp.size();
                }
                else
                {
                    // outer product [comp x output] -> temp
                    std::vector<double> temp(comp.size() * output_len, 0);
                    for (std::size_t j = 0; j < comp.size(); ++j)
                        for (std::size_t k = 0; k < output_len; ++k)
                            temp[j * output_len + k] = comp[j] * output[k];

                    // overwrite output with the newly calculated outer product
                    std::copy(temp.begin(), temp.end(), output);
                    output_len = temp.size();
                }
            }
        }

//        gauss_elim(n_terms, n_terms + 1, )
//        return ;
        return interp_matrix;
    }
};
