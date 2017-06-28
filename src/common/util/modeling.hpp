#ifndef SP2_MODELING_HPP
#define SP2_MODELING_HPP

#include "common/structure_t.hpp"

namespace sp2 {
namespace util {

/// create simple supercell by repeating input structure along it's lattice
/// vectors
/// \param input input primitive cell that is repeated along it's lattice
///        vectors to create a supercell
/// \param na number of times to repeat the structure in the direction of the a
///        lattice vector
/// \param nb number of times to repeat the structure in the direction of the b
///        lattice vector
/// \param nc number of times to repeat the structure in the direction of the c
///        lattice vector
/// \return constructed supercell
structure_t construct_supercell(const structure_t &input,
    int na, int nb = 1, int nc = 1);

structure_t graphene_unit_cell();

structure_t construct_arm_gnr(int width, int length, bool periodic = true);
structure_t construct_zz_gnr(int width, int length, bool periodic = true);

structure_t make_hydrogen_terminated(const structure_t &input);

structure_t center_by_avg(const structure_t &input);

structure_t make_capped_nanotube(const structure_t &cap,
    const structure_t &segment, int n_segment);

void make_nanotube_dataset();

} // namespace util
} // namespace sp2

#endif // SP2_MODELING_HPP
