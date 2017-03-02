#include "common/util/rotations.hpp"
#include <cmath>

using namespace std;
using namespace sp2;

void make_rotation_matrix(const vec3_t &axis, double cos_theta, double sin_theta,
    double output[3][3])
{
    // output = [u^T . u] * (1 - cos(theta))
    //        + [I] * cos(theta)
    //        + [u]_x * sin(theta)
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
            output[i][j] = (1 - cos_theta) * axis[i] * axis[j];

        output[i][i] += cos_theta;
        output[i][(i + 1) % 3] -= sin_theta * axis[(i + 2) % 3];
        output[i][(i + 2) % 3] += sin_theta * axis[(i + 1) % 3];
    }
}

void util::gen_rotation(vec3_t axis, double theta, double output[3][3])
{
    double cos_theta, sin_theta;
    sincos(theta, &sin_theta, &cos_theta);

    make_rotation_matrix(axis, cos_theta, sin_theta, output);
}


//void generate_orthobasis(double input[3][3])
//{
//    // Frisvad, Jeppe Revall. "Building an orthonormal basis from a 3D unit
//    // vector without normalization."
//    // Journal of Graphics Tools 16.3 (2012): 151-159.
//
//    const double &x = input[0][0],
//        &y = input[0][1],
//        &z = input[0][2];
//
//    // Handle the singularity
//    if (z < std::nextafter(-1.0, 0.0))
//    {
//        // matrix looks like
//        // input = {
//        //     { 0, 0,-1},
//        //     { 0,-1, 0},
//        //     {-1, 0, 0}
//        // };
//
//        input[1][1] = -1;
//        input[2][0] = -1;
//    }
//    else
//    {
//        double a = 1.0 / (1.0 + z),
//            b = -x * y * a;
//
//        // matrix looks like
//        // double output[3][3] = {
//        //     {              x,               y,  z},
//        //     {1.0 - x * x * a,               b, -x},
//        //     {              b, 1.0 - y * y * a, -y}
//        // };
//
//        input[1][0] = 1.0 - x * x * a,;
//        input[1][1] = b;
//        input[1][2] = -x;
//
//        input[2][0] = b;
//        input[2][1] = 1.0 - y * y * a;
//        input[2][2] = -y;
//    }
//}

//void util::gen_rotation(const vec3_t &a, const vec3_t &b, double output[3][3])
//{
////    auto axis = cross(a, b);
////
////    // a and b are unit vectors so we can do this
////    double sin_theta = axis.mag(),
////        cos_theta = dot(a, b);
////
////    // note, sin_theta is the magnitude of the cross product/axis
////    // so we just divide the axis vector by it to normalize it
////    make_rotation_matrix(axis / sin_theta, cos_theta, sin_theta, output);
//
//    double r1[3][3] = {
//        {a.x(), a.y(), a.z()},
//        {    0,     0,     0},
//        {    0,     0,     0}
//    };
//
//    generate_orthobasis(r1);
//
//    double r2[3][3] = {
//        {b.x(), b.y(), b.z()},
//        {    0,     0,     0},
//        {    0,     0,     0}
//    };
//
//    generate_orthobasis(r2);
//}

void util::gen_rand_rotation(double output[3][3])
{
    auto current_dir = vec3_t(0, 0, 1),
            new_dir  = random_vec3(),
            rot_axis = unit_normal_to(current_dir, new_dir);

    double rot1[3][3] = {},
            rot2[3][3] = {},
            theta1 = 2.0 * M_PI * rand() / (RAND_MAX + 1.0),
            theta2 = angle(current_dir, new_dir);

    // first rotate about z
    util::gen_rotation(current_dir, theta1, rot1);
    // then rotate to the randomly generated unit vector
    util::gen_rotation(rot_axis, theta2, rot2);

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            output[i][j] = 0;
            for (int k = 0; k < 3; ++k)
                output[i][j] += rot1[i][k] * rot2[k][j];
        }
    }
}