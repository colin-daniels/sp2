#include "common/minimize/minimize.hpp"

using namespace std;
using namespace sp2;

/// simple POD to hold information about linesearch boundaries
struct ls_bound_t
{
    double a, ///< alpha
        v, ///< value
        s; ///< slope
};

/// get the analytical minimum of a cubic fit to two input bounds (can return NaN)
double cubic_min(ls_bound_t low, ls_bound_t high);

/// get the analytical mininum of a quadratic fit to two input bounds
double quad_min(ls_bound_t low, ls_bound_t high);

double minimize::linesearch(oned_fn_t objective_fn, oned_fn_t slope_fn,
    double alpha)
{
    return linesearch([&](double a) {
        return std::make_pair(objective_fn(a), slope_fn(a));
    }, alpha);
}

double minimize::linesearch(diff1d_fn_t objective_fn, double alpha)
{
    constexpr int iter_lim = 16;

    // get values for alpha = 0
    double value, slope;
    tie(value, slope) = objective_fn(0);

    // record initial value
    auto initial_value = value;

    // right hand side quantities for the wolfe condition linesearch
    auto armijo = 1e-4 * slope, // sufficient decrease, from the armijo condition
        curvature = 1e-1 * slope; // the curvature condition

    // lower and upper bounds for minumum finding
    ls_bound_t low = {0, value, slope},
        high = {0, 0, 0};

    // running minimum, initialize with the information from alpha = 0
    auto min_point = make_pair(value, 0.0);

    for (int iteration = 0; iteration < iter_lim; ++iteration)
    {
        // check for errors in alpha
        if (!isfinite(alpha))
            return min_point.second;

        // update value and slope
        tie(value, slope) = objective_fn(alpha);

        // check the wolfe conditions
        if (value <= initial_value + alpha * armijo) // armijo
            if (abs(slope) <= abs(curvature))        // curvature
                return alpha;

        // update running minimum
        min_point = min(min_point, make_pair(value, alpha));

        // update the bounding interval for the minimum
        if (value < low.v && slope < 0 && alpha < high.a)
            low = {alpha, value, slope};
        else
            high = {alpha, value, slope};

        // get the new alpha
        auto minimum = cubic_min(low, high);
        if (isnormal(minimum))
            alpha = minimum;
        else
            alpha = quad_min(low, high);
    }

    // return the alpha that gave us the lowest value
    return min_point.second;
}

double cubic_min(ls_bound_t low, ls_bound_t high)
{
    auto delta_a = (high.a - low.a);

    // get cubic coefficients
    // f(x) = ax^3 + bx^2 + cx + d
    auto coeff_a = ((high.s + low.s) * delta_a - 2 * (high.v - low.v))
                   / (delta_a * delta_a * delta_a),
        coeff_b = (3 * (high.v - low.v) - (high.s + 2 * low.s) * delta_a)
                  / (delta_a * delta_a),
        coeff_c = low.s;

    // use the + root for the minimum position
    // (due to the positive second derivative)
    return low.a + (static_cast<double>(
                        sqrt(coeff_b * coeff_b - 3 * coeff_a * coeff_c))
                    - coeff_b) / (3 * coeff_a);
}

double quad_min(ls_bound_t low, ls_bound_t high)
{
    auto delta_a = (high.a - low.a);

    // get quadratic coefficients
    // f(x) = bx^2 + cx + d
    auto coeff_b = (high.v - low.v - low.s * delta_a) / (delta_a * delta_a),
        coeff_c = low.s;

    // return the calculated minimum
    return low.a - coeff_c / (2 * coeff_b);
}
