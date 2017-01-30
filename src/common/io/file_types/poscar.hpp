#ifndef SP2_POSCAR_HPP
#define SP2_POSCAR_HPP

#include <common/structure_t.hpp>

namespace sp2 {
namespace io {

bool read_poscar(std::string filename, structure_t &output);

bool write_poscar(std::string filename, const structure_t &input);

} // namespace io
} // namespace sp2

#endif // SP2_POSCAR_HPP
