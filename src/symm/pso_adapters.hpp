#ifndef SP2_PSO_ADAPTERS_HPP
#define SP2_PSO_ADAPTERS_HPP

#include "common/structure_t.hpp"
#include <vector>

namespace sp2 {
namespace symm {

/// places atoms normally, dim = n * 3 + 1
structure_t basic_adapter(
    const std::vector<double> &input,   ///< pso input
    double uc_range[2]                  ///< unit cell range (min/max)
);

std::vector<double> inverse_basic_adapter(
    const structure_t &structure,
    double uc_range[2]
);

/// places atoms as neighbors to other atoms, dim = n * 3 + 1
structure_t bonded_adapter(
    const std::vector<double> &input,   ///< pso input
    double uc_range[2],                 ///< unit cell range (min/max)
    double bond_range[2],               ///< atom bond range (min/max)
    std::size_t n_connected             ///< # of atoms placed independently
);

} // namespace symm
} // namespace sp2

#endif // SP2_PSO_ADAPTERS_HPP
