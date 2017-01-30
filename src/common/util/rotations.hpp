#ifndef SP2_MATH_HPP
#define SP2_MATH_HPP

#include "common/vec3_t.hpp"

namespace sp2 {
namespace util {

////////////////////////////////////////////////////////////////////////////////
// Miscellaneous Math Utility Functions                                       //
////////////////////////////////////////////////////////////////////////////////

/// generate a rotation matrix given an axis (unit vector) and angle
void gen_rotation(vec3_t axis, double theta, double output[3][3]);

/// generate a rotation matrix from the unit vector a to the unit vector b
void gen_rotation(const vec3_t &a, const vec3_t &b, double output[3][3]);

/// generate random rotation matrix
void gen_rand_rotation(double output[3][3]);

} // namespace util
} // namespace sp2

#endif // SP2_MATH_HPP
