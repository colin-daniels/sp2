#ifndef SP2_MODELING_HPP
#define SP2_MODELING_HPP

#include "common/structure_t.hpp"

namespace sp2 {
namespace util {

/// create simple supercell by repeating input structure along it's lattice
/// vectors
/// \param input input primitive cell that is repeated along it's lattice
///        vectors to create a supercell
/// \param supercell_dim supercell dimensions, aka how many times to repeat the
///        input structure in each periodic direction
/// \return constructed supercell
structure_t construct_supercell(const structure_t &input, int supercell_dim[3]);

structure_t construct_graphene(int m, int n, bool x_aligned = false);

structure_t make_hydrogen_terminated(const structure_t &input);

} // namespace util
} // namespace sp2

#endif // SP2_MODELING_HPP
