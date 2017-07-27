#ifndef SP2_STRUCTURAL_MUTATION_HPP
#define SP2_STRUCTURAL_MUTATION_HPP

// This belongs to metropolis.  However, apparently it is illegal to specialize
// sp2::enum_map from within sp2::minimize, so it has its own file to avoid the
// whole namespace song and dance.

#include "common/enums.hpp"
#include "common/python/numpy_util.hpp"

namespace sp2 {

/// types of changes that can be made to a structure
enum class structural_mutation_type : int
{
    INVALID = 0,
    CART_COORDS = 1,
    FRAC_COORDS = 2,
    LATTICE = 3
};

template<>
constexpr enum_map_t<structural_mutation_type>
    enum_map<structural_mutation_type> = {
    {structural_mutation_type::CART_COORDS, "carts"},
    {structural_mutation_type::FRAC_COORDS, "fracs"},
    {structural_mutation_type::LATTICE,     "lattice"}
};

/// result of a mutation function on a structure.
struct structural_mutation_t
{
    structural_mutation_type type;
    as_ndarray_t<double> data;
};

} // namespace sp2

#endif //SP2_STRUCTURAL_MUTATION_HPP
