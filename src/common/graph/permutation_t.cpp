#include "common/graph/permutation_t.hpp"

#include <queue>
#include <cstdlib>

using namespace std;
using namespace sp2;

graph::permutation_t::permutation_t(int n) :
    permut(n, 0)
{
    for (int i = 0; i < n; ++i)
        permut[i] = i;
}

graph::permutation_t::permutation_t(std::vector<int> permut_in) :
    permut(permut_in) {}

void graph::permutation_t::apply(ud_graph_t &graph) const
{
    vector<pair<int, int>> new_edges;
    new_edges.reserve(graph.n_edges() / 2);

    for (auto edge : graph.edges())
        if (edge.a < edge.b)
            new_edges.emplace_back(permut[edge.a], permut[edge.b]);

    graph = ud_graph_t(new_edges, graph.n_vertices());
}

graph::permutation_t graph::permutation_t::inverse() const
{
    vector<int> inv(permut.size());
    for (size_t i = 0; i < permut.size(); ++i)
        inv[permut[i]] = i;

    return permutation_t(inv);
}

void graph::permutation_t::randomize() {
    random_shuffle(permut.begin(), permut.end());}

int graph::get_bandwidth(const ud_graph_t &graph)
{
    int max_dist = 0;
    for (auto edge : graph.edges())
        max_dist = max(max_dist, static_cast<int>(std::abs(edge.a - edge.b)));
    return max_dist;
}

graph::permutation_t graph::cuthill_mckee(const ud_graph_t &graph)
{
    vector<int> offsets = graph.get_int_offsets(),
        edge_ids = graph.get_int_edge_ids();

    if (offsets.size() < 2)
        return permutation_t(offsets.size());

    int min_index = 0;
    vector<int> degree(offsets.size() - 1, 0),
        min_indices;

    for (size_t i = 0; i + 1 < offsets.size(); ++i)
    {
        degree[i] = offsets[i + 1] - offsets[i];
        if (i == 0 || degree[i] < degree[min_index])
        {
            min_index = i;
            min_indices = vector<int>(1, i);
        }
        else if (degree[i] == degree[min_index])
            min_indices.push_back(i);
    }

    min_index = min_indices[rand() % min_indices.size()];

    vector<int> permut,
        visited(offsets.size() - 1, 0);

    queue<int> to_visit;
    visited[min_index] = 1;
    to_visit.push(min_index);

    while (!to_visit.empty())
    {
        int init_size = to_visit.size();
        vector<pair<int, int>> visit_next;

        for (int i = 0; i < init_size; ++i)
        {
            int id = to_visit.front();
            permut.push_back(id);
            to_visit.pop();

            for (int j = offsets[id]; j < offsets[id + 1]; ++j)
            {
                int id_b = edge_ids[j];
                if (visited[id_b] == 0)
                {
                    visited[id_b] = 1;
                    visit_next.push_back(pair<int, int>(degree[id_b], id_b));
                }
            }
        }

        sort(visit_next.begin(), visit_next.end());
        for (size_t i = 0; i < visit_next.size(); ++i)
            to_visit.push(visit_next[i].second);
    }

    return graph::permutation_t(permut);
}

