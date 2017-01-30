#ifndef CGEN_EXP_HPP
#define CGEN_EXP_HPP

#include <vector>
#include <string>
#include <functional>

std::pair<double, double> max_error(const unsigned int n_points,
    const double low, const double high,
    std::function<double(double)> calc,
    std::function<double(double)> actual
);

/// generate lookup table for func() in the range [low:high)
std::vector<double> generate_lookup_table(unsigned int n_points,
    double low, double high, std::function<double(double)> func);

#endif // CGEN_EXP_HPP
