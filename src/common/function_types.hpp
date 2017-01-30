//
// Created by cc on 9/7/16.
//

#ifndef SP2_FUNCTION_TYPES_HPP
#define SP2_FUNCTION_TYPES_HPP

#include <functional>
#include <vector>

namespace sp2 {
////////////////////////////////////////////////////////////////////////////////
// Primary function types                                                     //
////////////////////////////////////////////////////////////////////////////////

/// R^n -> R^n (takes a vector, returns a vector, e.g. gradient)
using vector_fn_t = std::function<
        const std::vector<double>&(const std::vector<double>&)>;

/// R^n -> R (takes a vector, returns a double)
using scalar_fn_t = std::function<
        double(const std::vector<double>&)>;

/// R -> R (one dimensional function, takes a double, returns a double)
using oned_fn_t = std::function<double(double)>;

/// (R^n -> R, R^n -> R^n) Differentiable function, takes a vector and returns
/// a pair where the first element is the value and second is the gradient
using diff_fn_t = std::function<
        std::pair<double, std::vector<double>>(const std::vector<double>&)>;

/// (R -> R, R -> R) One dimensional differentiable function, takes a vector and
/// returns a pair where the first element is the value and second is the slope
using diff1d_fn_t = std::function<std::pair<double, double>(double)>;

/// void function
template<class T>
using void_fn_t = std::function<void(T)>;

////////////////////////////////////////////////////////////////////////////////
// Utilities                                                                  //
////////////////////////////////////////////////////////////////////////////////

oned_fn_t to_one_dim(
    scalar_fn_t input_fn,
     const std::vector<double> &initial_pos,
     const std::vector<double> &direction
);

oned_fn_t to_one_dim(
    scalar_fn_t input_fn,
    const std::vector<double> &direction
);

diff1d_fn_t to_one_dim(
    diff_fn_t input_fn,
   const std::vector<double> &initial_pos,
   const std::vector<double> &direction
);

diff1d_fn_t to_one_dim(
    diff_fn_t input_fn,
    const std::vector<double> &direction
);

} // namespace sp2

#endif //SP2_FUNCTION_TYPES_HPP
