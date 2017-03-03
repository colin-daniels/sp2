#ifndef SP2_UTIL_MATH_HPP
#define SP2_UTIL_MATH_HPP

#include <cinttypes>

namespace sp2 {

/// constexpr power function for integer powers
template<class T = std::uint64_t>
constexpr T int_pow(T base, std::uint64_t exp)
{
    std::int64_t result = 1;
    while(exp)
    {
        if(exp & 1)
            result *= base;

        exp >>= 1;
        base *= base;
    }

    return result;
}

/// constexpr factorial
constexpr std::uint64_t factorial(std::uint64_t n)
{
    std::uint64_t result = 1;
    for (std::uint64_t i = 1; i <= n; ++i)
        result *= i;

    return result;
}

/// get the factorial "difference", high! / low!
constexpr std::uint64_t factorial_diff(std::uint64_t high, std::uint64_t low)
{
    std::uint64_t result = 1;
    for (std::uint64_t i = low + 1; i <= high; ++i)
        result *= i;

    return result;
}

} // namespace sp2

#endif // SP2_UTIL_MATH_HPP
