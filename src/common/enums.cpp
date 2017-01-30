#define BASIC_ENUMS_HPP_IMPL
#include "common/enums.hpp"
#undef BASIC_ENUMS_HPP_IMPL

#include "io/util.hpp"
#include <unordered_map>

using namespace std;

namespace sp2 {

// atom types
template<> std::string enum_to_str(atom_type type) {
    return type == atom_type::HYDROGEN ? "H" : "C";}

template<> atom_type enum_from_str(std::string str) {
    return str == "H" ? atom_type::HYDROGEN : atom_type::CARBON;}

} // namespace scpp