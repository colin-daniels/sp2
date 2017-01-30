#include "outputs/outputs.hpp"
#include "util/stream_io.hpp"
#include "util/lookup_table.hpp"

#include "common/util/math.hpp"
#include "common/util/timing.hpp"
#include "common/util/random.hpp"

#include <iostream>
#include <vector>
#include <cmath>

namespace {
/// calculate the exponential lookup table
std::vector<double> calc_exp_table(int n_points, double range)
{
    double step = range / n_points;

    std::vector<double> result(n_points, 0);
    for (int i = 0; i < n_points; ++i)
        result[i] = std::exp((i + 0.5) * step);

    return result;
}

constexpr int n_lk = 8192;
constexpr double range_lk = -4.7204523127 * 2.1;
const std::vector<double> exp_table = calc_exp_table(n_lk, range_lk);
} // anonymous namespace

/// exponential lookup table function, valid in the range [-9.9:0]
__attribute__((always_inline)) static inline double lk_exp(double x)
{
    constexpr auto step = range_lk / n_lk,
        inv_step = n_lk / range_lk;

    const double dst = x - static_cast<int>(inv_step * x) * step - step / 2;
    return (1.0 + dst + (dst * dst) / 2)
           * exp_table[static_cast<int>(inv_step * x)];
}

void exp_table_hpp()
{
    std::ofstream outfile("exp_table.hpp");
    outfile.precision(std::numeric_limits<double>::max_digits10);

    outfile << "constexpr double exp_table_step = "
            << (range_lk / n_lk) << ";\n"
            << "constexpr double exp_table[] = ";
    write_array(outfile, exp_table);
    outfile << ";" << std::endl;
    outfile.close();
}

#ifdef SP2_DEBUG
#include <gtest/gtest.h>

TEST(codegen, benchmark) {

    sp2::util::rng_t rng;
    std::vector<double> data(1000);
    for (auto &d : data)
        d = rng.rand(range_lk, 0.0);

    unsigned int n_iters = 10'000;

    double lk_time = sp2::benchmark(n_iters, [&]{
        volatile double x = 0;
        for (auto d : data)
            x = lk_exp(d);
    });

    double exp_time = sp2::benchmark(n_iters, [&]{
        volatile double x = 0;
        for (auto d : data)
            x = std::exp(d);
    });

    auto lk_err = max_error(1'000'000, 0, range_lk,
        [](double x){
            return lk_exp(x);
        }, [](double x){
            return std::exp(x);
    });

    std::cout << "  lk time: " << lk_time  << '\n'
              << " exp time: " << exp_time << '\n'
              << "  speedup: " << (exp_time / lk_time) << '\n'
              << "max error: " << lk_err.first << " at x = " << lk_err.second
              << std::endl;
}

#endif // SP2_DEBUG