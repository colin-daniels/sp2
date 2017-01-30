#ifndef SP2_STRUCTURE_HPP
#define SP2_STRUCTURE_HPP

#include <string>
#include "common/structure_t.hpp"

namespace sp2 {

/// file types (mostly for structure output)
enum class file_type : int
{
    UNKNOWN = 0, ///< invalid type
    AUTO    = 1, ///< type should be inferred from filename
    XYZ     = 2,
    JSON    = 3,
    CIF     = 4,
    POSCAR  = 5
};

namespace io {
/// infer filetype from filename
file_type infer_filetype(std::string filename);

/// write structure to a file, determines format from filename
bool write_structure(std::string filename,
    const structure_t &input, bool append = false,
    file_type type = file_type::AUTO);

/// read structure from a file, determines format from filename
bool read_structure(std::string filename,
    structure_t &output,
    file_type type = file_type::AUTO);

/// convert one (structure) file type to another, data loss possible
bool convert_file(std::string input_file, std::string output_file,
    file_type input_type = file_type::AUTO,
    file_type output_type = file_type::AUTO);

} // namespace io
} // namespace sp2


#endif // SP2_STRUCTURE_HPP
