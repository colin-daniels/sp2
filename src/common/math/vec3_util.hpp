#ifndef SP2_VEC3_UTIL_HPP
#define SP2_VEC3_UTIL_HPP

#include "common/math/vec3_t.hpp"

#include <vector>
#include <limits>

namespace sp2 {

/// vector dot product
constexpr double dot(const vec3_t &a, const vec3_t &b);
/// vector cross product
constexpr vec3_t cross(const vec3_t &a, const vec3_t &b);
/// calculate the angle between two vectors (uses std::atan2)
inline double angle(const vec3_t &a, const vec3_t &b);

/// generate a randomly oriented unit vector using the provided
/// UniformRandomBitGenerator

/// generate a randomly oriented unit vector using a provided random generator
/// \tparam URBG generator type, satisfies UniformRandomBitGenerator
/// \param g random bit generator of type URBG
/// \return randomly oriented unit vector
template<class URBG>
vec3_t random_vec3(URBG &&g);
/// generate a randomly oriented unit vector using std::rand()
inline vec3_t random_vec3();

/// generate a unit vector normal to a single input vector
inline vec3_t unit_normal_to(const vec3_t &a);
/// generate a unit vector normal to two input vectors
inline vec3_t unit_normal_to(const vec3_t &a, const vec3_t &b);

/// get minimum elements of two vectors
inline vec3_t min_elem(const vec3_t &a, const vec3_t &b);
/// get maximum elements of two vectors
inline vec3_t max_elem(const vec3_t &a, const vec3_t &b);

/// get the bounds of a set of vec3_t's (aka min/max x, y, z)
inline std::pair<vec3_t, vec3_t> get_bounds(const std::vector<vec3_t> &input);

/// convert vector of vec3_t into doubles (x1, y1, z1, x2, y2, z2, ...)
inline std::vector<double> v3tod(const std::vector<vec3_t> &input);
/// convert vector of doubles (x1, y1, z1, x2, y2, z2, ...) to vec3_t
inline std::vector<vec3_t> dtov3(const std::vector<double> &input);

////////////////////////////////////////////////////////////////////////////////
// Implementations                                                            //
////////////////////////////////////////////////////////////////////////////////

constexpr double dot(const vec3_t &a, const vec3_t &b)
{
    return a.x * b.x
         + a.y * b.y
         + a.z * b.z;
}

constexpr vec3_t cross(const vec3_t &a, const vec3_t &b)
{
    return vec3_t(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    );
}

inline double angle(const vec3_t &a, const vec3_t &b)
{
    return std::atan2(cross(a, b).mag(), dot(a, b));
}

inline vec3_t unit_normal_to(const vec3_t &a)
{
    // need +1 on z to make sure the second vector isn't 0 and isn't equal to a
    return unit_normal_to(a, {a.x, a.y, a.z + 1});
}

inline  vec3_t unit_normal_to(const vec3_t &a, const vec3_t &b)
{
    return cross(a, b).unit_vector();
}

inline vec3_t min_elem(const vec3_t &a, const vec3_t &b)
{
    return vec3_t(
        std::min(a.x, b.x),
        std::min(a.y, b.y),
        std::min(a.z, b.z)
    );
}

inline vec3_t max_elem(const vec3_t &a, const vec3_t &b)
{
    return vec3_t(
        std::max(a.x, b.x),
        std::max(a.y, b.y),
        std::max(a.z, b.z)
    );
}

inline std::pair<vec3_t, vec3_t> get_bounds(const std::vector<vec3_t> &input)
{
    using nl = std::numeric_limits<double>;
    auto result = std::make_pair(
        vec3_t{   nl::max(),    nl::max(),    nl::max()}, // min
        vec3_t{nl::lowest(), nl::lowest(), nl::lowest()}  // max
    );

    for (const auto& v : input)
        result = {min_elem(result.first, v),
                  max_elem(result.second, v)};

    return result;
}

template<class URBG>
vec3_t random_vec3(URBG &&g)
{
    constexpr auto rng_min = std::remove_reference_t<URBG>::min(),
        rng_max = std::remove_reference_t<URBG>::max();

    constexpr double inv_range = 1.0 /
        static_cast<double>(rng_max - rng_min);

    const double theta = 2 * M_PI * inv_range * (g() - rng_min), // [0:2pi]
        u = 1 - 2 * inv_range * (g() - rng_min); // [-1:1]

    return vec3_t{u,
        std::sqrt(1 - u * u) * std::cos(theta),
        std::sqrt(1 - u * u) * std::sin(theta)
    };
}

inline vec3_t random_vec3()
{
    struct rand_struct
    {
        using result_type = decltype(std::rand());

        static constexpr result_type max() { return RAND_MAX; }
        static constexpr result_type min() { return 0; }
        result_type operator()() { return std::rand(); }
    };

    return random_vec3(
        rand_struct{}
    );
}

inline std::vector<double> v3tod(const std::vector<vec3_t> &input)
{
    std::vector<double> output;
    output.reserve(input.size() * 3);
    for (auto &v : input)
        for (auto &d : v)
            output.push_back(d);

    return output;
}

inline std::vector<vec3_t> dtov3(const std::vector<double> &input)
{
    std::vector<vec3_t> output;
    output.reserve(input.size() / 3);
    for (size_t i = 0; i < input.size(); i += 3)
        output.emplace_back(&input[i]);

    return output;
}

} // namespace sp2

#endif // SP2_VEC3_UTIL_HPP
