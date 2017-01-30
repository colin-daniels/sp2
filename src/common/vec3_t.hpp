#ifndef SP2_VEC3_T_HPP
#define SP2_VEC3_T_HPP

#include <cmath>
#include <utility>
#include <vector>

#include <boost/mpi/datatype_fwd.hpp>
#include <boost/mpl/bool.hpp>

// for boost::mpi interaction with vec3_t
namespace boost {
namespace serialization {
class access;
} // namespace serialization
} // namespace boost


namespace sp2 {
namespace util {
class rng_t;
} // namespace util

////////////////////////////////////////////////////////////////////////////////
// Utility R^3 vector class                                                   //
////////////////////////////////////////////////////////////////////////////////
struct alignas(16) vec3_t
{
    double pos[3];

    vec3_t() = default;

    explicit vec3_t(const double *ptr) noexcept :
            pos{ptr[0], ptr[1], ptr[2]}
    {}

    vec3_t(double x_in, double y_in, double z_in) noexcept :
            pos{x_in, y_in, z_in}
    {}

    operator double *()
    { return pos; }

    operator const double *() const
    { return pos; }

////////////////////////////////////////////////////////////////////////////////
// Access-related functions                                                   //
////////////////////////////////////////////////////////////////////////////////

    double *begin()
    { return &pos[0]; }

    double *end()
    { return &pos[3]; }

    const double *begin() const
    { return &pos[0]; }

    const double *end() const
    { return &pos[3]; }

    constexpr std::size_t size() const
    { return 3; }

    double x() const
    { return pos[0]; }

    double y() const
    { return pos[1]; }

    double z() const
    { return pos[2]; }

    double &operator[](std::size_t i)
    { return pos[i]; }

    const double &operator[](std::size_t i) const
    { return pos[i]; }

////////////////////////////////////////////////////////////////////////////////
// Mathematical operators                                                     //
////////////////////////////////////////////////////////////////////////////////

    vec3_t operator+(const vec3_t &o) const
    {
        return vec3_t(x() + o.x(), y() + o.y(), z() + o.z());
    }

    vec3_t &operator+=(const vec3_t &other)
    {
        *this = *this + other;
        return *this;
    }

    vec3_t operator-(const vec3_t &o) const
    {
        return vec3_t(x() - o.x(), y() - o.y(), z() - o.z());
    }

    vec3_t operator-() const
    {
        return -1 * *this;
    }

    vec3_t &operator-=(const vec3_t &other)
    {
        *this = *this - other;
        return *this;
    }

    friend vec3_t operator*(const vec3_t &v, const double a)
    {
        return vec3_t(v.x() * a, v.y() * a, v.z() * a);
    }

    friend vec3_t operator*(const double a, const vec3_t &v)
    {
        return v * a;
    }

    vec3_t &operator*=(const double a)
    {
        *this = *this * a;
        return *this;
    }

    friend vec3_t operator/(const vec3_t &v, const double a)
    {
        return vec3_t(v.x() / a, v.y() / a, v.z() / a);
    }

    vec3_t &operator/=(const double a)
    {
        *this = *this / a;
        return *this;
    }

////////////////////////////////////////////////////////////////////////////////
// Miscellaneous member functions                                             //
////////////////////////////////////////////////////////////////////////////////

    /// returns the magnitude squared for the vector
    double mag_sq() const
    {
        return x() * x() + y() * y() + z() * z();
    }

    /// returns the 2-norm of the vector
    double mag() const
    {
        return static_cast<double>(std::sqrt(mag_sq()));
    }

    /// set equal to a random unit vector (use rand())
    vec3_t &randomize();

    /// set equal to a random unit vector (use provided random number generator)
    vec3_t &randomize(util::rng_t &rng);

    /// normalize the vector
    vec3_t &normalize()
    {
        *this /= mag();
        return *this;
    }

    /// return the vectors unit vector does not modify it
    vec3_t unit_vector() const
    { return *this / mag(); }

    /// multiply a 3x3 matrix by this vector
    vec3_t mul_3x3(const double mat[3][3]) const
    {
        return vec3_t(
                x() * mat[0][0] + y() * mat[0][1] + z() * mat[0][2],
                x() * mat[1][0] + y() * mat[1][1] + z() * mat[1][2],
                x() * mat[2][0] + y() * mat[2][1] + z() * mat[2][2]
        );
    }

    // boost::serialization (for boost::mpi)
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive &ar, const unsigned int) const
    {
        ar & pos;
    }
};

static_assert(std::is_pod<vec3_t>::value, "");

/// vector dot product
inline double dot(const vec3_t &a, const vec3_t &b)
{
    return a.x() * b.x()
           + a.y() * b.y()
           + a.z() * b.z();
}

/// vector cross product
inline vec3_t cross(const vec3_t &a, const vec3_t &b)
{
    return vec3_t(
            a.y() * b.z() - a.z() * b.y(),
            a.z() * b.x() - a.x() * b.z(),
            a.x() * b.y() - a.y() * b.x()
    );
}

/// calculate the angle between two vectors
inline double angle(const vec3_t &a, const vec3_t &b)
{
    return static_cast<double>(std::atan2(cross(a, b).mag(), dot(a, b)));
}

/// generate a vector normal to a single input vector
vec3_t unit_normal(const vec3_t &a);

/// generate a vector normal to both of the input vectors
inline vec3_t unit_normal(const vec3_t &a, const vec3_t &b)
{
    return cross(a, b).unit_vector();
}

/// get minimum elements of two vectors
inline vec3_t min_elem(const vec3_t &a, const vec3_t &b)
{
    return vec3_t(
        std::min(a[0], b[0]),
        std::min(a[1], b[1]),
        std::min(a[2], b[2])
    );
}

/// get maximum elements of two vectors
inline vec3_t max_elem(const vec3_t &a, const vec3_t &b)
{
    return vec3_t(
        std::max(a[0], b[0]),
        std::max(a[1], b[1]),
        std::max(a[2], b[2])
    );
}

// TODO: fix swap
//constexpr void swap(vec3_t &a, vec3_t &b)
//{
//    vec3_t temp = a;
//    a = b;
//    b = a;
//}

////////////////////////////////////////////////////////////////////////////////
// Utility Functions                                                          //
////////////////////////////////////////////////////////////////////////////////

/// convert vector of vec3_t into doubles (x1, y1, z1, x2, y2, z2, ...)
std::vector<double> v3tod(const std::vector<vec3_t> &input);

/// convert vector of doubles (x1, y1, z1, x2, y2, z2, ...) to vec3_t
std::vector<vec3_t> dtov3(const std::vector<double> &input);

} // namespace sp2

// for boost::mpi
namespace boost {
namespace mpi {

template<>
struct is_mpi_datatype<sp2::vec3_t> :
        public mpl::true_ {};

} // namespace mpi
} // namespace boost

#endif
