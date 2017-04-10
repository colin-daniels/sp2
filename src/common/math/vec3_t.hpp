#ifndef SP2_VEC3_T_HPP
#define SP2_VEC3_T_HPP

#include <cmath>
#include <cstdlib>
#include <cstddef>
#include <utility>

namespace sp2 {

/// used in SFINAE to determine which types can participate in division
/// or multiplication of vec3_t objects
template<class T>
constexpr bool vec3_math_compat_v = std::is_arithmetic<std::decay_t<T>>::value;

////////////////////////////////////////////////////////////////////////////////
// Utility R^3 vector class                                                   //
////////////////////////////////////////////////////////////////////////////////
struct alignas(16) vec3_t
{
//    double x, y, z;
    union
    {
        struct {
            double x, y, z;
        };
        double data[3];
    };

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

    constexpr double *begin() {return data;}
    constexpr double *end() {return begin() + 3;}

    constexpr const double *begin() const {
        return data;}
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

    constexpr vec3_t& operator+=(const vec3_t &o) {
        x += o.x; y += o.y; z += o.z; return *this; }

    constexpr vec3_t& operator-=(const vec3_t &o) {
        x -= o.x; y -= o.y; z -= o.z; return *this; }

    template<class T, class = std::enable_if_t<vec3_math_compat_v<T>>>
    constexpr vec3_t& operator*=(T a) {
        x *= a; y *= a; z *= a; return *this; }

    template<class T, class = std::enable_if_t<vec3_math_compat_v<T>>>
    constexpr vec3_t& operator/=(T a) {
        x /= a; y /= a; z /= a; return *this; }


    constexpr vec3_t operator-() const {
        return -1.0 * *this;}

    constexpr vec3_t operator+(const vec3_t &o) const {
        return {x + o.x, y + o.y, z + o.z};}

    constexpr vec3_t operator-(const vec3_t &o) const {
        return {x - o.x, y - o.y, z - o.z};}


    template<class T, class = std::enable_if_t<vec3_math_compat_v<T>>>
    constexpr vec3_t operator*(T a) const { return {x * a, y * a, z * a}; }

    template<class T, class = std::enable_if_t<vec3_math_compat_v<T>>>
    constexpr vec3_t operator/(T a) const { return {x / a, y / a, z / a}; }


    template<class T, class = std::enable_if_t<vec3_math_compat_v<T>>>
    friend constexpr vec3_t operator*(T&& a, const vec3_t &v) {
        return v * std::forward<T>(a); }

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

// TODO: fix
//static_assert(vec3_t{1, 2, 3}[0] == 1, "");
//static_assert(vec3_t{1, 2, 3}[1] == 2, "");
//static_assert(vec3_t{1, 2, 3}[2] == 3, "");

} // namespace sp2

#ifdef SP2_ENABLE_MPI

#include <boost/mpi/datatype_fwd.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/serialization/is_bitwise_serializable.hpp>

// for boost::mpi
namespace boost {

namespace serialization {
// vec3_t is a trivial type so it can be serialized bitwise
template<>
struct is_bitwise_serializable<sp2::vec3_t> :
    mpl::true_ {};
} // namespace serialization

namespace mpi {
// again because vec3_t is trivial, it is also an mpi datatype
template<>
struct is_mpi_datatype<sp2::vec3_t> :
    mpl::true_ {};
} // namespace mpi

} // namespace boost

#endif // SP2_ENABLE_MPI

#endif // SP2_VEC3_T_HPP
