#ifndef SP2_SYMM_SYSTEM_CONTROL_T_HPP
#define SP2_SYMM_SYSTEM_CONTROL_T_HPP


#include "common/structure_t.hpp"
#include "common/math/vec3_t.hpp"
#include "airebo/system_control_t.hpp"
#include "common/neighbor/bond_control_t.hpp"
#include "symm/space_group_t.hpp"

#include <vector>

namespace sp2 {
namespace symm {

class system_control_t
{
private:
    double unit_cell;
    space_group_t group;

    double total_potential;

    double lattice[3][3] = {};

    std::vector<atom_type_old> types;
    std::vector<vec3_t> position,
        gradient;

    fbc::bond_control_t bond_control;
    airebo::system_control_t sub_sys;

public:
    system_control_t();

    void init(const space_group_t &space_group,
        const structure_t &structure);

    /// set system space group (and appropriate symmetries)
    void set_group(const space_group_t &input);
    /// get system structure (irreducible representation)
    void set_structure(const structure_t &input);

    /// set system structure (irreducible representation)
    structure_t get_structure() const;
    /// get the system structure (full structure with symmetries applied)
    structure_t get_full_structure() const;

    /// get gradient
    std::vector<double> get_gradient() const;
    /// get total potential
    double get_value() const;

    /// calculate forces and total potential for the current structure
    void update();

private:
    /// get the 'full' structure but limit ourselves to atoms within max_dist
    /// bonds of the atoms in the irreducible representation
    structure_t get_limited_structure(size_t max_dist = 3);

    std::vector<vec3_t> get_full_pos() const;
    std::vector<atom_type_old> get_full_types() const;
};

} // namespace symm
} // namespace sp2

#endif // SP2_SYMM_SYSTEM_CONTROL_T_HPP
