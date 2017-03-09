#ifndef SP2_PHONOPY_IO_HPP
#define SP2_PHONOPY_IO_HPP

#include <string>
#include <vector>
#include <common/structure_t.hpp>
#include <common/math/vec3_t.hpp>

namespace sp2 {
namespace phonopy {

std::vector<std::string> read_irreps(std::string filename = "irreps.yaml");

std::tuple<structure_t, std::vector<std::pair<int, vec3_t>>>
    read_displacements(std::string filename = "disp.yaml");

void write_force_sets(std::string filename);

} // namespace phonopy
} // namespace sp2

#endif // SP2_PHONOPY_IO_HPP
