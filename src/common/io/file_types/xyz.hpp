#ifndef SP2_IO_XYZ_HPP
#define SP2_IO_XYZ_HPP

/// \file xyz.hpp
/// \brief Header for .XYZ file IO functions

#include <string>
#include <vector>

namespace sp2 {
namespace io {

/// read a frame of an xyz file
bool read_xyz(std::string filename,
    std::string &comment_line,
    std::vector<std::string> &types,
    std::vector<double> &positions,
    int frame = 0);

/// write an xyz file (or append to one)
bool write_xyz(std::string filename,
    const std::string &comment_line,
    const std::vector<std::string> &types,
    const std::vector<double> &positions,
    bool append = false);

/// \brief get the number of frames in an xyz file
/// \returns -1 if an error is encountered
int get_xyz_frames(std::string filename);

} // namespace io
} // namespace sp2

#endif // SP2_IO_XYZ_HPP
