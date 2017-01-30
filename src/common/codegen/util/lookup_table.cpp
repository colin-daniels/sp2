#include "lookup_table.hpp"
#include <limits>
#include <iostream>
#include <utility>
#include <cmath>

std::pair<double, double> max_error(const unsigned int n_points,
    const double low, const double high,
    std::function<double(double)> calc,
    std::function<double(double)> actual
)
{
    const double step = (high - low) / n_points;

    auto max_err = std::make_pair(0.0, low);
    for (auto i = 0U; i < n_points; ++i)
    {
        auto x = i * step + low;
        auto actual_val = actual(x);
        max_err = std::max(max_err, std::make_pair(
            std::abs((calc(x) - actual_val) / actual_val), x));
    }

    return max_err;
}

std::vector<double> generate_lookup_table(unsigned int n_points,
    double low, double high, std::function<double(double)> func)
{
    const unsigned int test_points = 100;
    const double range = high - low,
        step = range / n_points;

    std::vector<double> values(n_points + 1);

    auto error = [&](double x) {
        return std::abs(values[static_cast<int>((x - low) / step)] - func(x));
    };

    auto get_error = [&](double t_low, double t_step) {
        double max_error = 0.0;
        for (auto i = 0U; i < test_points; ++i)
            max_error = std::max(max_error, error(i * t_step + t_low));

        return max_error;
    };

    double max_error = 0;
    for (auto i = 0U; i < n_points; ++i)
    {
        double min_error = std::numeric_limits<double>::max(),
            min_value = 0;

        double t_low = i * step + low,
            t_step = step / test_points;

        for (auto j = 0U; j < test_points; ++j)
        {
            values[i] = func(t_step * j + t_low);
            auto err = get_error(t_low, t_step);

            if (err < min_error)
            {
                min_value = values[i];
                min_error = err; //
            }
        }

        values[i] = min_value;
        max_error = std::max(max_error, min_error);
    }

    std::cout << "lookup table max error: " << max_error << std::endl;
    return values;
}
