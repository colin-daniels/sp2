#ifndef SP2_IO_CIF_HPP
#define SP2_IO_CIF_HPP

#include <string>
#include "common/structure_t.hpp"

namespace sp2 {
namespace io {

bool write_cif(std::string filename, const structure_t &structure);

} // namespace io
} // namespace sp2

#endif // SP2_IO_CIF_HPP
