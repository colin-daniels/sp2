#ifndef SP2_MATH_HPP
#define SP2_MATH_HPP

#include "common/math/vec3_t.hpp"
#include "common/math/mat3x3_t.hpp"

namespace sp2 {
namespace util {

////////////////////////////////////////////////////////////////////////////////
// Miscellaneous Math Utility Functions                                       //
////////////////////////////////////////////////////////////////////////////////

/// generate a rotation matrix given an axis (unit vector) and angle
mat3x3_t gen_rotation(vec3_t axis, double theta);

/// generate a rotation matrix from the unit vector a to the unit vector b
mat3x3_t gen_rotation(const vec3_t &from, const vec3_t &to);

} // namespace util
} // namespace sp2

#endif // SP2_MATH_HPP
