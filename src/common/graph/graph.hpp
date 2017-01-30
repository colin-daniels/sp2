#ifndef SP2_GRAPH_HPP
#define SP2_GRAPH_HPP

/// \file graph.hpp
/// \brief Utility functions involving graphs, particularly undirected unweighted ones

#include <vector>
#include <string>

#include "common/graph/ud_graph_t.hpp"

namespace sp2 {
/// graph function namespace
namespace graph {

/// breadth first search from start->end, returns path from end->start
/// with length less than max_dist (no limit if max_dist <= 0)
std::vector<int> bfs_path(const ud_graph_t &graph,
    int start, int end, int max_dist);

/// breadth first search from start->end, returns distance less than
/// max_dist (no limit if max_dist <= 0) or -1 if no path is found
int bfs_dist(const ud_graph_t &graph,
    int start, int end, int max_dist);

/// depth first search from start->end, returns path from end->start
/// with length less than max_dist (no limit if max_dist <= 0)
std::vector<int> dfs_path(const ud_graph_t &graph,
    int start, int end, int max_dist);

/// depth first search from start->end, returns distance less than
/// max_dist (no limit if max_dist <= 0) or -1 if no path is found
int dfs_dist(const ud_graph_t &graph,
    int start, int end, int max_dist);

/// get rings
/// \param max_ring int maximum size of rings to look for, omit or -1 for no limit
/// \return vector of rings (vectors)
std::vector<std::vector<int>> bfs_rings(const ud_graph_t &graph,
    int max_ring = -1);

/// write adjacency matrix to file
void write_gplot_graph(std::string filename, const ud_graph_t &graph);

} // namespace graph
} // namespace sp2

#endif // SP2_GRAPH_HPP
