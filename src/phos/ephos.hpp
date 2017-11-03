#pragma once

#include "src/common/function_types.hpp"
#include "src/common/structure_t.hpp"
#include "src/common/math/vec3_t.hpp"
#include "src/common/neighbor/bond_control_t.hpp"
#include "src/common/graph/ud_graph_t.hpp"

#include <vector>

namespace sp2 {
namespace phos {

class phosphorene_sys_t
{
public:
    explicit phosphorene_sys_t(const sp2::structure_t &input);

    void update();

    /// get the system potential
    double get_value() const;

    structure_t get_structure() const { return structure; }

    /// get atom positions
    std::vector<double> get_position() const;

    /// get atom forces
    std::vector<double> get_gradient() const;

    /// set atom positions
    /// \param input const vector<double>& input atom positions
    void set_position(const std::vector<double> &input);

    diff_fn_t get_diff_fn();

private:
    static constexpr double bond_cutoff = 3.0; // Angstroms

    structure_t structure;

    double potential;
    std::vector<vec3_t> forces;

    sp2::fbc::bond_control_t bond_control;
    std::vector<vec3_t> pristine_deltas;
};

} // namespace phos
} // namespace sp2
