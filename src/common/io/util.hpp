#ifndef SP2_UTIL_HPP
#define SP2_UTIL_HPP

#include <string>
#include <vector>

namespace sp2 {
namespace io {

/// clear the contents of a file
bool clear_file(std::string filename);
/// copy a file from source to dest
bool copy_file(std::string source, std::string dest);
/// check the existence of a file
bool file_exists(std::string filename);
/// remove a file
bool remove_file(std::string filename);
/// read a file directly into a string
bool read_file(std::string filename, std::string &content);
/// write a string directly to a file
bool write_file(std::string filename, const std::string &content);

/// trim string
void trim(std::string &str, std::string chars = " ");
/// split line by characters, resulting strings do not contain any of the
/// characters that were specified to be split on
std::vector<std::string> split_line(const std::string &line,
    const std::string &characters = " \n\r\t");

/// convert string to uppercase
std::string toupper(const std::string &input);
/// convert string to lowercase
std::string tolower(const std::string &input);

} // namespace io
} // namespace sp2

#endif // SP2_UTIL_HPP
