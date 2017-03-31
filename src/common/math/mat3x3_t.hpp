#ifndef SP2_MAT3X3_T_HPP
#define SP2_MAT3X3_T_HPP

#include <utility>

namespace sp2 {

struct mat3x3_t
{
public:
    double data[3][3];

    mat3x3_t() = default;
    constexpr mat3x3_t(const double (&data)[3][3]) :
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

    constexpr       double *begin()       { return std::begin(data[0]); }
    constexpr const double *begin() const { return std::begin(data[0]); }

    constexpr       double *end()       { return std::end(data[2]); }
    constexpr const double *end() const { return std::end(data[2]); }

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

    /// non-modifying inverse, returns a new matrix
    constexpr mat3x3_t inverse() const
    {
        // (1 / det) *
        //
        // |22 21| |02 01| |12 11|
        // |12 11| |22 21| |02 01|
        //
        // |20 22| |00 02| |10 12|
        // |10 12| |20 22| |00 02|
        //
        // |21 20| |01 00| |11 10|
        // |11 10| |21 20| |01 00|
        //
        auto& d = data;
        mat3x3_t inverse_mat{{
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

        auto det = inverse_mat[0][0] * inverse_mat[0][0] +
                   inverse_mat[1][1] * inverse_mat[1][1] +
                   inverse_mat[2][2] * inverse_mat[2][2];

        for (auto &v : inverse_mat)
            v /= det;

        return inverse_mat;
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
        return vec.mul_3x3(data);
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
