#ifndef SP2_SYMM_UTIL_HPP
#define SP2_SYMM_UTIL_HPP

#include <vector>
#include <map>
#include <string>

#include "symm/space_group_t.hpp"
#include "common/structure_t.hpp"

namespace sp2 {
namespace symm {

void test();

/// force directing atoms towards the gyroid surface
double gyroid_ptnl(const double *pos, double *deriv, const double c = 0);

// structure_t get_irreducible(const structure_t &input,
//     const std::map<std::string, space_group_t> &groups);

} // namespace symm
} // namespace sp2

#endif // SP2_SYMM_UTIL_HPP
