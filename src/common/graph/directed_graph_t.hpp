#ifndef SP2_DIRECTED_GRAPH_T_HPP
#define SP2_DIRECTED_GRAPH_T_HPP

#include <vector>
#include <unordered_set>

namespace sp2 {
namespace graph {

/// directed graph type, designed for fast construction and access, but
/// modification is relatively slow
///
/// \tparam Edge
template<class Edge>
class directed_graph_t
{
private:
    std::vector<Edge> edge_vector;

    std::vector<std::size_t> vertex_indices;
public:
    directed_graph_t() = default;

    decltype(auto) edges() const { return edge_vector; }
    decltype(auto) edges(std::size_t vertex_idx) const
    {
        return make_range(
            edge_vector.begin() + vertex_indices[vertex_idx],
            edge_vector.begin() + vertex_indices[vertex_idx + 1]
        );
    }

    decltype(auto) begin() const { return edge_vector.begin(); }
    decltype(auto) end()   const { return edge_vector.end(); }

    std::size_t n_edges()    const { return edge_vector.size(); }
    std::size_t n_vertices() const { return vertex_indices.size() - 1; }

};

} // namespace graph
} // namespace sp2

#endif // SP2_DIRECTED_GRAPH_T_HPP
