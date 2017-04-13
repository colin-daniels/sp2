#ifndef SP2_NUMERICAL_DIFF_HPP
#define SP2_NUMERICAL_DIFF_HPP

#include <functional>
#include <array>

namespace sp2 {
namespace util {

template<std::size_t N, class T>
using cd_coeff_t = std::pair<T, std::array<T, N>>;

/// calculate numerical derivative (1D) via central difference
template<std::size_t N = 9, class T = double>
constexpr auto central_difference(std::function<T(T)> func, T x_0, T step_size)
{
    constexpr auto coeff = std::get<cd_coeff_t<N, T>>(
        std::make_tuple(
            // format is {denominator coefficient, {numerator coefficients}}
            cd_coeff_t<3, T>(  3, {{1, 0, -1}}),
            cd_coeff_t<5, T>( 12, {{1, -8, 0, 8, -1}}),
            cd_coeff_t<7, T>( 60, {{-1, 9, -45, 0, 45, -9, 1}}),
            cd_coeff_t<9, T>(840, {{3, -32, 168, -672, 0, 672, -168, 32, -3}})
        )
    );

    auto sum = T{},
            x = x_0 - step_size * (N / 2);
    for (std::size_t i = 0; i < N; ++i, x += step_size)
        if (i != N / 2)
            sum += coeff.second[i] * func(x);

    return sum / (coeff.first * step_size);
}

} // namespace util
} // namespace sp2

#endif // SP2_NUMERICAL_DIFF_HPP
