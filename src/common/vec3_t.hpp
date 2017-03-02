#ifndef SP2_VEC3_T_HPP
#define SP2_VEC3_T_HPP

#include <cmath>
#include <cstdlib>
#include <cstddef>
#include <utility>
#include <vector>

namespace sp2 {
namespace util {
class rng_t;
} // namespace util

////////////////////////////////////////////////////////////////////////////////
// Utility R^3 vector class                                                   //
////////////////////////////////////////////////////////////////////////////////
struct vec3_t;

/// vector dot product
constexpr double dot(const vec3_t &a, const vec3_t &b);
/// vector cross product
constexpr vec3_t cross(const vec3_t &a, const vec3_t &b);
/// calculate the angle between two vectors (uses std::atan2)
inline double angle(const vec3_t &a, const vec3_t &b);

/// generate a randomly oriented unit vector using the provided generator
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

/// Utility R^3 vector class
struct alignas(16) vec3_t
{
    double x, y, z;

    vec3_t() = default;

    constexpr explicit vec3_t(const double *ptr) :
        x(ptr[0]), y(ptr[1]), z(ptr[2]) {}

    constexpr vec3_t(double x, double y, double z) :
        x(x), y(y), z(z) {}

    constexpr operator const double*() const { return begin(); }
    constexpr operator       double*()       { return begin(); }

////////////////////////////////////////////////////////////////////////////////
// Access-related functions                                                   //
////////////////////////////////////////////////////////////////////////////////

    constexpr double *begin() {return reinterpret_cast<double*>(this);}
    constexpr double *end() {return begin() + 3;}

    constexpr const double *begin() const {
        return reinterpret_cast<const double*>(this);}
    constexpr const double *end() const {
        return begin() + 3;}

    constexpr double& operator[](std::size_t i) {
        return begin()[i];}
    constexpr const double& operator[](std::size_t i) const {
        return begin()[i];}

    constexpr std::size_t size() const {return 3;}

////////////////////////////////////////////////////////////////////////////////
// Mathematical operators                                                     //
////////////////////////////////////////////////////////////////////////////////

    constexpr vec3_t operator-() const {
        return -1.0 * *this;}

    constexpr vec3_t operator+(const vec3_t &o) const {
        return {x + o.x, y + o.y, z + o.z};}

    constexpr vec3_t operator-(const vec3_t &o) const {
        return {x - o.x, y - o.y, z - o.z};}


    template<class T>
    friend constexpr vec3_t operator*(T&& a, const vec3_t &v)
    {
        return {
            std::forward<T>(a) * v.x,
            std::forward<T>(a) * v.y,
            std::forward<T>(a) * v.z
        };
    }

    template<class T>
    friend constexpr vec3_t operator*(const vec3_t &v, T&& a) {
        return std::forward<T>(a) * v;}

    template<class T>
    friend constexpr vec3_t operator/(const vec3_t &v, T&& a)
    {
        return {
            v.x / std::forward<T>(a),
            v.y / std::forward<T>(a),
            v.z / std::forward<T>(a)
        };
    }


    constexpr vec3_t& operator+=(const vec3_t &other) {
        return *this = *this + other;}

    constexpr vec3_t& operator-=(const vec3_t &other) {
        return *this = *this - other;}

    template<class T>
    constexpr vec3_t& operator*=(T&& a) {
        return *this = *this * std::forward<T>(a);}

    template<class T>
    constexpr vec3_t& operator/=(T&& a) {
        return *this = *this / std::forward<T>(a);}


////////////////////////////////////////////////////////////////////////////////
// Miscellaneous member functions                                             //
////////////////////////////////////////////////////////////////////////////////

    /// returns the magnitude squared for the vector
    constexpr double mag_sq() const {
        return x * x + y * y + z * z;}

    /// returns the 2-norm of the vector
    double mag() const {
        return std::sqrt(mag_sq());}

    /// normalize the vector
    vec3_t &normalize() {
        return *this /= mag();}

    /// return the vectors unit vector does not modify it
    vec3_t unit_vector() const {
        return *this / mag();}

    /// multiply a 3x3 matrix by this vector
    constexpr vec3_t mul_3x3(const double (&mat)[3][3]) const
    {
        vec3_t result{0, 0, 0};
        return vec3_t{
            x * mat[0][0] + y * mat[0][1] + z * mat[0][2],
            x * mat[1][0] + y * mat[1][1] + z * mat[1][2],
            x * mat[2][0] + y * mat[2][1] + z * mat[2][2]
        };
    }
};

static_assert(std::is_trivial<vec3_t>::value, "");
static_assert(sizeof(vec3_t) == sizeof(double[4]), "");

static_assert(offsetof(vec3_t, x) == 0, "");
static_assert(offsetof(vec3_t, y) == alignof(double), "");
static_assert(offsetof(vec3_t, z) == 2 * alignof(double), "");

static_assert(vec3_t{1, 2, 3}.x == 1, "");
static_assert(vec3_t{1, 2, 3}.y == 2, "");
static_assert(vec3_t{1, 2, 3}.z == 3, "");

static_assert(vec3_t{1, 2, 3}[0] == 1, "");
static_assert(vec3_t{1, 2, 3}[1] == 2, "");
static_assert(vec3_t{1, 2, 3}[2] == 3, "");

/// vector dot product
constexpr double dot(const vec3_t &a, const vec3_t &b)
{
    return a.x * b.x
         + a.y * b.y
         + a.z * b.z;
}

/// vector cross product
constexpr vec3_t cross(const vec3_t &a, const vec3_t &b)
{
    return vec3_t(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    );
}

/// calculate the angle between two vectors
inline double angle(const vec3_t &a, const vec3_t &b)
{
    return std::atan2(cross(a, b).mag(), dot(a, b));
}

/// generate a unit vector normal to a single input vector
inline vec3_t unit_normal_to(const vec3_t &a)
{
    // need +1 on z to make sure the second vector isn't 0 and isn't equal to a
    return unit_normal_to(a, {a.x, a.y, a.z + 1});
}

/// generate a unit vector normal to both of the input vectors
inline  vec3_t unit_normal_to(const vec3_t &a, const vec3_t &b)
{
    return cross(a, b).unit_vector();
}

/// get minimum elements of two vectors
inline vec3_t min_elem(const vec3_t &a, const vec3_t &b)
{
    return vec3_t(
        std::min(a.x, b.x),
        std::min(a.y, b.y),
        std::min(a.z, b.z)
    );
}

/// get maximum elements of two vectors
inline vec3_t max_elem(const vec3_t &a, const vec3_t &b)
{
    return vec3_t(
        std::max(a.x, b.x),
        std::max(a.y, b.y),
        std::max(a.z, b.z)
    );
}

template<class URBG>
vec3_t random_vec3(URBG &&g)
{
    constexpr double inv_range = 1.0 /
        static_cast<double>(URBG::max() - URBG::min());

    const double theta = 2 * M_PI * inv_range * (g() - URBG::min()), // [0:2pi]
        u = 1 - 2 * inv_range * (g() - URBG::min()); // [-1:1]

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

////////////////////////////////////////////////////////////////////////////////
// Utility Functions                                                          //
////////////////////////////////////////////////////////////////////////////////

/// convert vector of vec3_t into doubles (x1, y1, z1, x2, y2, z2, ...)
std::vector<double> v3tod(const std::vector<vec3_t> &input);

/// convert vector of doubles (x1, y1, z1, x2, y2, z2, ...) to vec3_t
std::vector<vec3_t> dtov3(const std::vector<double> &input);

} // namespace sp2

#include <boost/mpi/datatype_fwd.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/serialization/is_bitwise_serializable.hpp>

// for boost::mpi
namespace boost {

namespace serialization {
template<>
struct is_bitwise_serializable<sp2::vec3_t> :
    mpl::true_ {};
} // namespace serialization

namespace mpi {
template<>
struct is_mpi_datatype<sp2::vec3_t> :
    mpl::true_ {};
} // namespace mpi

} // namespace boost

#endif
