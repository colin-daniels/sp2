#ifndef SP2_SYMM_SETTINGS_T_HPP
#define SP2_SYMM_SETTINGS_T_HPP

#include "common/json/json_serializable_t.hpp"
#include "common/minimize/settings.hpp"

#include <string>

namespace sp2 {
namespace symm {

/// structural information, atom positions/types/lattice vectors
struct symm_settings_t : public io::json_serializable_t
{
    /// number of atoms (irreducible representation)
    int n_atoms = 0;

    /// space group name
    std::string space_group = "I a -3 d",
    /// file from which to read the space groups
        space_group_file = "space-groups.txt";

    /// min/max unit cell size for the search
    double unit_cell_range[2] = {},
    /// min/max distance to place atoms
        bond_range[2] = {1.1, 1.7};

    /// max number of disconnected subgraphs (in terms of bond structure)
    int n_connected = 1;

    /// where to output the structures found
    std::string output_filename = "output.xyz";

    /// particle swarm settings
    minimize::pso_settings_t pso_set;

    /// whether to use conjugate gradient in addition to the pso
    bool use_cg = true;
    /// conjugate gradient settings
    minimize::acgsd_settings_t acgsd_set;

    symm_settings_t()
    {
        acgsd_set.output_level = 0;
        acgsd_set.iteration_limit = 100;
    }

    bool serialize(Json::Value &output) const;
    bool deserialize(const Json::Value &input);
};

} // namespace symm
} // namespace sp2

#endif // SP2_SYMM_SETTINGS_T_HPP
