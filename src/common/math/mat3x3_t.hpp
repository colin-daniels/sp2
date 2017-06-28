#ifndef SP2_MAT3X3_T_HPP
#define SP2_MAT3X3_T_HPP

#include <utility>
#include "vec3_t.hpp"
#include "vec3_util.hpp"

namespace sp2 {

class mat3x3_t
{
private:
    vec3_t data[3];

public:
    mat3x3_t() = default;
    constexpr mat3x3_t(const vec3_t &a, const vec3_t &b, const vec3_t &c) :
        data{a, b, c} {}

    explicit constexpr mat3x3_t(const double (&data)[3][3]) :
        data{
            {data[0][0], data[0][1], data[0][2]},
            {data[1][0], data[1][1], data[1][2]},
            {data[2][0], data[2][1], data[2][2]}
        } {}

    constexpr operator const decltype(data)&() const { return data; }
    constexpr operator       decltype(data)&()       { return data; }

////////////////////////////////////////////////////////////////////////////////
// Access-related functions                                                   //
////////////////////////////////////////////////////////////////////////////////

    constexpr       auto begin()       { return std::begin(data); }
    constexpr const auto begin() const { return std::begin(data); }

    constexpr       auto end()       { return std::end(data); }
    constexpr const auto end() const { return std::end(data); }

    constexpr decltype(auto) operator[](int i) const { return data[i]; }
    constexpr decltype(auto) operator[](int i)       { return data[i]; }

    /// modifying transpose, transposes current matrix and returns a reference
    constexpr mat3x3_t& transpose()
    {
        for (int i = 0; i < 3; ++i)
        {
            for (int j = i + 1; j < 3; ++j)
            {
                double tmp = data[i][j];
                data[i][j] = data[j][i];
                data[j][i] = tmp;
            }
        }

        return *this;
    }

    /// non-modifying transpose, returns a new matrix
    constexpr mat3x3_t transposed() const
    {
        return mat3x3_t{{
            {data[0][0], data[1][0], data[2][0]},
            {data[0][1], data[1][1], data[2][1]},
            {data[0][2], data[1][2], data[2][2]}
        }};
    }

    // This is 'inline' to decrease the likelihood of missed
    // opportunities for dead code elimination in determinant().
    /// get the adjugate matrix
    constexpr inline mat3x3_t adjugate() const
    {
        auto& d = data;
        return mat3x3_t{{
            {(d[2][2] * d[1][1]) - (d[1][2] * d[2][1]),
             (d[0][2] * d[2][1]) - (d[2][2] * d[0][1]),
             (d[1][2] * d[0][1]) - (d[0][2] * d[1][1])},

            {(d[2][0] * d[1][2]) - (d[1][0] * d[2][2]),
             (d[0][0] * d[2][2]) - (d[2][0] * d[0][2]),
             (d[1][0] * d[0][2]) - (d[0][0] * d[1][2])},

            {(d[2][1] * d[1][0]) - (d[1][1] * d[2][0]),
             (d[0][1] * d[2][0]) - (d[2][1] * d[0][0]),
             (d[1][1] * d[0][0]) - (d[0][1] * d[1][0])}
        }};
    }

    /// non-modifying inverse, returns a new matrix
    constexpr mat3x3_t inverse() const
    {
        auto adj = adjugate();
        auto cofacs = adj.transposed();
        double det = dot(cofacs.data[0], data[0]);

        for (auto &v : adj.data)
            v /= det;
        return adj;
    }

    /// compute the determinant
    constexpr double determinant() const
    {
        // For the determinant, take any single row of the matrix, and dot
        // it with the corresponding row from the cofactor matrix.

        // This implementation relies on dead code elimination for the two
        // thirds of adjugate() that we don't use.
        //
        // optimized x86_64 assembly looks approximately like this:
        //
        //     https://gist.github.com/ExpHP/3d24e243ad65e2c7cfcddc6116f03d93
        //
        auto cofacs = adjugate().transposed();
        return dot(cofacs.data[0], data[0]);
    }

    /// modifying inverse, inverts the current matrix and returns a reference
    constexpr mat3x3_t& invert()
    {
        return *this = inverse();
    }

    /// get the 3x3 identity matrix
    static constexpr mat3x3_t identity()
    {
        return mat3x3_t{{
            {1, 0, 0},
            {0, 1, 0},
            {0, 0, 1}
        }};
    }

    sp2::vec3_t operator*(const sp2::vec3_t &vec) const
    {
        return {
            dot(data[0], vec),
            dot(data[1], vec),
            dot(data[2], vec)
        };
    }

    sp2::mat3x3_t operator*(const sp2::mat3x3_t &m) const
    {
        sp2::mat3x3_t result;
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                result[i][j] = 0;
                for (int k = 0; k < 3; ++k)
                    result[i][j] += data[i][k] * m[k][j];
            }
        }

        return result;
    }

    sp2::mat3x3_t& operator+=(const sp2::mat3x3_t &other)
    {
        for (int i = 0; i < 3; ++i)
            data[i] += other.data[i];

        return *this;
    }
};

static_assert(std::is_trivial<mat3x3_t>::value, "");

} // namespace sp2


#ifdef SP2_ENABLE_MPI

#include <boost/mpi/datatype_fwd.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/serialization/is_bitwise_serializable.hpp>

// similar to sp2::vec3_t, we need to tell boost that mat3x3_t is trivial
// as well as an mpi datatype
namespace boost {

namespace serialization {
template<>
struct is_bitwise_serializable<sp2::mat3x3_t> :
    mpl::true_ {};
} // namespace serialization

namespace mpi {
template<>
struct is_mpi_datatype<sp2::mat3x3_t> :
    mpl::true_ {};
} // namespace mpi

} // namespace boost

#endif // SP2_ENABLE_MPI

#endif // SP2_MAT3X3_T_HPP
