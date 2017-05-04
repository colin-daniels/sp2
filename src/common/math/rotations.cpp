#include "rotations.hpp"
#include "common/math/vec3_util.hpp"
#include "mat3x3_t.hpp"
#include <cmath>

using namespace std;
using namespace sp2;

sp2::mat3x3_t make_rotation_matrix(const vec3_t &axis,
    double cos_theta, double sin_theta)
{
    // output = [u^T . u] * (1 - cos(theta))
    //        + [I] * cos(theta)
    //        + [u]_x * sin(theta)
    sp2::mat3x3_t rotation;
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
            rotation[i][j] = (1 - cos_theta) * axis[i] * axis[j];

        rotation[i][i] += cos_theta;
        rotation[i][(i + 1) % 3] -= sin_theta * axis[(i + 2) % 3];
        rotation[i][(i + 2) % 3] += sin_theta * axis[(i + 1) % 3];
    }

    return rotation;
}

mat3x3_t util::gen_rotation(vec3_t axis, double theta)
{
    double cos_theta, sin_theta;
    sincos(theta, &sin_theta, &cos_theta);

    return make_rotation_matrix(axis, cos_theta, sin_theta);
}

constexpr sp2::mat3x3_t generate_orthobasis(const sp2::vec3_t &input)
{
    // Frisvad, Jeppe Revall. "Building an orthonormal basis from a 3D unit
    // vector without normalization."
    // Journal of Graphics Tools 16.3 (2012): 151-159.

    // Handle the singularity (note: "close enough" to -1.0, more than
    // one ulp away though)
    if (input.z() < -0.99999999999)
    {
        // output basis is
        return {
            { 0,  0, -1}, // input
            { 0, -1,  0},
            {-1,  0,  0}
        };
    }
    else
    {
        const double a = 1.0 / (1.0 + input.z()),
            b = -input.x() * input.y() * a;

        // output basis is
        return {
            {input.x(), input.y(), input.z()}, // input
            {/*x*/ 1.0 - input.x() * input.x() * a, /*y*/ b, /*z*/ -input.x()},
            {/*x*/ b, /*y*/ 1.0 - input.y() * input.y() * a, /*z*/ -input.y()}
        };
    }
}

// expects a.mag() = 1, b.mag() = 1
mat3x3_t util::gen_rotation(const vec3_t &a, const vec3_t &b)
{
    // make orthogonal bases for a and b
    sp2::mat3x3_t basis_a = generate_orthobasis(a),
        basis_b = generate_orthobasis(b);

    // R = A*B^T now, doing R * v yields a coordinate change from a to b for v
    mat3x3_t rotation;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            rotation[i][j] = dot(basis_a[i], basis_b[j]);

    return rotation;
}
