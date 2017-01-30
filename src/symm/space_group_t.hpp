#ifndef SP2_SPACE_GROUP_T_HPP
#define SP2_SPACE_GROUP_T_HPP

#include <array>
#include <vector>
#include <string>
#include <fstream>
#include <map>

#include "common/vec3_t.hpp"

namespace sp2 {
namespace symm {

/// simple class to store and apply symmetries corresponding to a space group
struct space_group_t
{
public:
    /// space group number
    int number = 0;
    /// space group name (H-M/intl notation)
    std::string name;
    /// actual transformation matrices (4x4)
    std::vector<std::array<double, 16>> transform;

    /// get the space group number
    int get_number() const;
    /// get the number of symmetries
    std::size_t n_symm() const;
    /// get the space group name (intl notation)
    std::string get_name() const;

    /// apply the space group's symmetries to an input vector
    /// and output the result
    std::vector<vec3_t> apply_symm(const std::vector<vec3_t> &input) const;

    /// parse from an input stream
    std::istream& parse(std::istream &stream);
};

/// Read space groups into a map from a specified file.
std::map<std::string, space_group_t> read_groups(std::string filename);

} // namespace symm
} // namespace sp2

#endif // SP2_SPACE_GROUP_T_HPP
