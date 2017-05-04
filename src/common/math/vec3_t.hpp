#ifndef SP2_VEC3_T_HPP
#define SP2_VEC3_T_HPP

#include <cmath>
#include <cstdlib>
#include <cstddef>
#include <iterator>
#include <utility>

namespace sp2 {

/// used in SFINAE to determine which types can participate in division
/// or multiplication of vec3_t objects
template<class T>
constexpr bool vec3_math_compat_v = std::is_arithmetic<std::decay_t<T>>::value;

////////////////////////////////////////////////////////////////////////////////
// Utility R^3 vector class                                                   //
////////////////////////////////////////////////////////////////////////////////
class alignas(32) vec3_t
{
private:
    double data[3];

public:

    vec3_t() = default;

    constexpr explicit vec3_t(const double *ptr) noexcept :
        data{ptr[0], ptr[1], ptr[2]} {}

    constexpr vec3_t(double x, double y, double z) noexcept :
        data{x, y, z} {}

////////////////////////////////////////////////////////////////////////////////
// Comparison operators                                                       //
////////////////////////////////////////////////////////////////////////////////

    friend constexpr bool operator==(const vec3_t &a, const vec3_t &b) noexcept
    {
        // note: bitwise and
        return (a[0] == b[0]) &
               (a[1] == b[1]) &
               (a[2] == b[2]);
    }

    friend constexpr bool operator!=(const vec3_t &a, const vec3_t &b) noexcept
    {
        return !(a == b);
    }

    // swap?

////////////////////////////////////////////////////////////////////////////////
// Access-related functions                                                   //
////////////////////////////////////////////////////////////////////////////////

    constexpr       double *begin()       {return std::begin(data);}
    constexpr const double *begin() const {return std::begin(data);}

    constexpr       double *end()       {return std::end(data);}
    constexpr const double *end() const {return std::end(data);}

    constexpr double& operator[](std::size_t i)       {return data[i];}
    constexpr double  operator[](std::size_t i) const {return data[i];}

    constexpr std::size_t size() const {return 3;}

    constexpr double& x()       {return data[0];}
    constexpr double  x() const {return data[0];}

    constexpr double& y()       {return data[1];}
    constexpr double  y() const {return data[1];}

    constexpr double& z()       {return data[2];}
    constexpr double  z() const {return data[2];}

////////////////////////////////////////////////////////////////////////////////
// Mathematical operators                                                     //
////////////////////////////////////////////////////////////////////////////////

    constexpr vec3_t& operator+=(const vec3_t &o)
    {
        for (int i = 0; i < 3; ++i)
            data[i] += o.data[i];
        return *this;
    }

    constexpr vec3_t& operator-=(const vec3_t &o)
    {
        for (int i = 0; i < 3; ++i)
            data[i] -= o.data[i];
        return *this;
    }

    template<class T, class = std::enable_if_t<vec3_math_compat_v<T>>>
    constexpr vec3_t& operator*=(T a)
    {
        for (auto& v : data)
            v *= a;
        return *this;
    }

    template<class T, class = std::enable_if_t<vec3_math_compat_v<T>>>
    constexpr vec3_t& operator/=(T a)
    {
        for (auto& v : data)
            v /= a;
        return *this;
    }

    constexpr vec3_t operator-() const { return *this * -1.0; }

    constexpr vec3_t operator+(const vec3_t &o) const
    {
        return {
            data[0] + o[0],
            data[1] + o[1],
            data[2] + o[2]
        };
    }

    constexpr vec3_t operator-(const vec3_t &o) const
    {
        return {
            data[0] - o[0],
            data[1] - o[1],
            data[2] - o[2]
        };
    }

    template<class T, class = std::enable_if_t<vec3_math_compat_v<T>>>
    constexpr vec3_t operator*(T a) const
    {
        return {
            data[0] * a,
            data[1] * a,
            data[2] * a
        };
    }

    template<class T, class = std::enable_if_t<vec3_math_compat_v<T>>>
    constexpr vec3_t operator/(T a) const
    {
        return {
            data[0] / a,
            data[1] / a,
            data[2] / a
        };
    }

    template<class T, class = std::enable_if_t<vec3_math_compat_v<T>>>
    friend constexpr vec3_t operator*(T&& a, const vec3_t &v) {
        return v * std::forward<T>(a); }


    /// multiply a 3x3 matrix by this vector (as a column vector)
    template<class T, class = std::enable_if_t<vec3_math_compat_v<T>>>
    friend constexpr vec3_t operator*(const T (&mat)[3][3], const vec3_t &v)
    {
        return {
            v[0] * mat[0][0] + v[1] * mat[0][1] + v[2] * mat[0][2],
            v[0] * mat[1][0] + v[1] * mat[1][1] + v[2] * mat[1][2],
            v[0] * mat[2][0] + v[1] * mat[2][1] + v[2] * mat[2][2]
        };
    }

////////////////////////////////////////////////////////////////////////////////
// Miscellaneous member functions                                             //
////////////////////////////////////////////////////////////////////////////////

    /// returns the magnitude squared for the vector
    constexpr double mag_sq() const {
        return data[0] * data[0]
             + data[1] * data[1]
             + data[2] * data[2];
    }

    /// returns the 2-norm of the vector
    double mag() const {
        return std::sqrt(mag_sq());}

    /// normalize the vector, assumes magnitude is not 0
    vec3_t &normalize() {
        return *this /= mag();}

    /// return the vectors unit vector, non-modifying
    vec3_t unit_vector() const {
        return *this / mag();}
};

static_assert(std::is_trivial<vec3_t>::value, "");
static_assert(sizeof(vec3_t) == 32, "");

static_assert(vec3_t{1, 2, 3}.x() == 1, "");
static_assert(vec3_t{1, 2, 3}.y() == 2, "");
static_assert(vec3_t{1, 2, 3}.z() == 3, "");

static_assert(vec3_t{1, 2, 3}[0] == 1, "");
static_assert(vec3_t{1, 2, 3}[1] == 2, "");
static_assert(vec3_t{1, 2, 3}[2] == 3, "");

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
