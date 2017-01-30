#ifndef SP2_PERMUTATION_T_HPP
#define SP2_PERMUTATION_T_HPP

#include <vector>

#include "common/graph/ud_graph_t.hpp"

namespace sp2 {
namespace graph {

class permutation_t
{
private:
    std::vector<int> permut; ///< permut[old_id] = new_id

public:
    explicit permutation_t(int n);
    permutation_t(std::vector<int> permut_in);

    permutation_t inverse() const;
    void apply(ud_graph_t &graph) const;

    void randomize();
};

/// Get the input graph's bandwidth
int get_bandwidth(const ud_graph_t &graph);

/// Cuthill-Mckee bandwidth reduction
permutation_t cuthill_mckee(const ud_graph_t &graph);

} // namespace graph
} // namespace sp2

#endif // SP2_PERMUTATION_T_HPP
