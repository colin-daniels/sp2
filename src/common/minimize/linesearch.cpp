#include "common/minimize/minimize.hpp"
#include "settings.hpp"

#include <utility>
#include <tuple>
#include <cmath>
#include <iostream>

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

class uphill_linesearch: public exception
{
    virtual const char* what() const noexcept
    { return "Uphill linesearch"; }
};

double minimize::linesearch(const ls_settings_t &settings,
    double alpha, oned_fn_t objective_fn, oned_fn_t slope_fn)
{
    return linesearch(settings, alpha, [&](double a) {
        return std::make_pair(objective_fn(a), slope_fn(a));
    });
}

double minimize::linesearch(const ls_settings_t &settings,
    double alpha, diff1d_fn_t objective_fn)
{
    // get values for alpha = 0
    double value, slope;
    std::tie(value, slope) = objective_fn(0);

    if (alpha <= 0)
        throw logic_error("linesearch initial alpha <= 0");
    if (slope > 0)
        throw uphill_linesearch();

    // TODO
    if (settings.weak_force_threshold != 0)
        throw runtime_error("weak force linesearch not yet implemented");

    auto initial_value = value;

    // Right hand side quantities for the wolfe condition linesearch.
    // - sufficient decrease
    auto armijo = settings.armijo_threshold * std::abs(slope);
    // - the curvature condition
    auto curvature = settings.curvature_threshold * std::abs(slope);

    // lower and upper bounds for minimum finding
    // (hard lower bound, soft upper bound)
    ls_bound_t low = {0, value, slope},
        high = {0, 0, 0};

    // running minimum, initialize with the information from alpha = 0
    auto min_point = std::make_pair(value, 0.0);

    for (int iteration = 0; iteration < settings.iteration_limit; ++iteration)
    {
        // check for errors in alpha
        if (!std::isfinite(alpha))
            return min_point.second;

        // update value and slope
        std::tie(value, slope) = objective_fn(alpha);

        // check the wolfe conditions
        if (value <= initial_value - alpha * armijo)    // armijo
            if (std::abs(slope) <= curvature) // curvature
                return alpha;

        // update running minimum
        min_point = std::min(min_point, std::make_pair(value, alpha));

        // update the bounding interval for the minimum
        if (value < low.v && slope < 0 && alpha < high.a)
            low = {alpha, value, slope};
        else
            high = {alpha, value, slope};

        // get the new alpha
        auto minimum = cubic_min(low, high);
        if (std::isnormal(minimum))
            alpha = minimum;
        else
            alpha = quad_min(low, high);
    }

    // return the alpha that gave us the lowest value
    // note: This could be zero!
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

    // f'(x) has two roots.
    // For real a, b, c, the + root always corresponds to the minimum.
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
