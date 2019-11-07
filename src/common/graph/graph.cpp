#include "common/graph/graph.hpp"
#include "common/minimize/minimize.hpp"

#include <boost/serialization/utility.hpp>

#include <queue>
#include <cmath>
#include <fstream>

#include <thread>
#include <future>
#include <numeric>

using namespace std;
using namespace sp2;

std::vector<int> graph::bfs_path(const ud_graph_t &graph,
    int start, int end, int max_dist)
{
    queue<int> to_visit;
    vector<int> visited(graph.n_vertices(), 0);

    to_visit.push(start);
    visited[start] = 1;

    int depth = 1;
    while (!to_visit.empty() && !visited[end])
    {
        int queue_size = to_visit.size();
        for (int i = 0; i < queue_size; ++i)
        {
            int vertex = to_visit.front();
            to_visit.pop();

            for (auto edge : graph.edges(vertex))
            {
                if (visited[edge.b])
                    continue;

                visited[edge.b] = depth + 1;
                to_visit.push(edge.b);
            }
        }

        if (max_dist > 0 && depth >= max_dist)
            break;
        else
            depth += 1;
    }

    if (!visited[end])
        return vector<int>();

    vector<int> path(1, end);
    for (int next_depth = depth; next_depth > 0; --next_depth)
    {
        int vertex = path.back();
        for (auto edge : graph.edges(vertex))
        {
            if (visited[edge.b] == next_depth)
            {
                path.push_back(edge.b);
                break;
            }
        }
    }

    return path;
}

int graph::bfs_dist(const ud_graph_t &graph,
    int start, int end, int max_dist)
{
    auto path = bfs_path(graph, start, end, max_dist);
    return static_cast<int>(path.size()) - 1;
}

std::vector<int> graph::dfs_path(const ud_graph_t &graph,
    int start, int end, int max_dist)
{
    vector<int> visited(graph.n_vertices(), false);

    vector<int> path;
    path.push_back(start);

    while (!path.empty() && path.back() != end)
    {
        int id = path.back();
        visited[id] = true;

        if (max_dist <= 0 || static_cast<int>(path.size()) < max_dist)
        {
            auto neigh = graph.neighbors(id);

            // find next vertex that has not been visited
            auto next = find_if(neigh.begin(), neigh.end(),
                [&visited](auto neighbor_id) {
                    return !visited[neighbor_id];
                });

            // if one was found we continue, and don't remove
            // anything from the current path vector
            if (next != neigh.end())
            {
                path.push_back(*next);
                continue;
            }
        }

        // if no new vertex was found, we move backwards
        path.pop_back();
    }

    return path;
}

int graph::dfs_dist(const ud_graph_t &graph,
    int start, int end, int max_dist)
{
    auto path = dfs_path(graph, start, end, max_dist);
    return static_cast<int>(path.size()) - 1;
}

void bfs_all_dist_wk(int root, const std::vector<int> &offsets, const std::vector<int> &edges, int max_dist, int *visited)
{
    queue<int> to_visit;
    to_visit.push(root);
    visited[root] = -1;

    for (int dist = 1; !to_visit.empty() && (max_dist == -1 || dist <= max_dist); ++dist)
    {
        int current_size = to_visit.size();
        for (int i = 0; i < current_size; ++i)
        {
            int id = to_visit.front();
            to_visit.pop();

            for (int j = offsets[id]; j < offsets[id + 1]; ++j)
            {
                int id_b = edges[j];
                if (visited[id_b] == 0)
                {
                    visited[id_b] = dist;
                    to_visit.push(id_b);
                }
            }
        }
    }
    visited[root] = 0;
}

std::vector<int> bfs_all_pairs(const std::vector<int> &offsets, const std::vector<int> &edges, int max_dist)
{
    int nv = offsets.size() - 1;
    vector<int> result(nv * nv, 0);
    for (int i = 0; i < nv; ++i)
        bfs_all_dist_wk(i, offsets, edges, max_dist, &result[i * nv]);

    return result;
}

void bfs_ring_check_path(vector<int> &path_a, vector<int> &path_b, vector<int> &all_pairs,
    const vector<int> &offsets, const vector<int> &edges)
{
    int nv = static_cast<int>(sqrt(all_pairs.size())),
        id_a = path_a.back(),
        id_b = path_b.back(),
        dist = all_pairs[id_a * nv + id_b];

    // exit condition
    if (static_cast<int>(path_a.size() + path_b.size()) > 2 * dist)
    {
        int rem_n = 0;
        if (path_a.back() == path_b.front())
        {
            path_a.pop_back();
            rem_n |= 1;
        }

        if (path_b.back() == path_a.front())
        {
            path_b.pop_back();
            rem_n |= 2;
        }

        // make sure rings actually form a loop
        if (rem_n == 3)
            return;
        else if (rem_n > 0)
        {
            int end_a, end_b;
            if (rem_n == 1)
            {
                end_a = path_a.front();
                end_b = path_b.back();
            }
            else
            {
                end_a = path_a.back();
                end_b = path_b.front();
            }
            for (int i = offsets[end_a]; i < offsets[end_a + 1]; ++i)
                if (edges[i] == end_b)
                    return;
        }

        path_a.clear();
        path_b.clear();

        return;
    }

    int last_a = path_a.size() > 1 ? *(path_a.end() - 2) : -1,
        last_b = path_b.size() > 1 ? *(path_b.end() - 2) : -1;

    // find next pair in the ring
    for (int i = offsets[id_a]; i < offsets[id_a + 1]; ++i)
    {
        int id_c = edges[i];

        // dont go backwards
        if (id_c == last_a)
            continue;

        for (int j = offsets[id_b]; j < offsets[id_b + 1]; ++j)
        {
            int id_d = edges[j];
            if (id_d == last_b)
                continue;

            int cd_dist = all_pairs[id_c * nv + id_d],
                ad_dist = all_pairs[id_a * nv + id_d],
                bc_dist = all_pairs[id_b * nv + id_c];

            // next pair needs to be the same distance
            if (cd_dist != dist || ad_dist > dist || bc_dist > dist ||
                ad_dist == 0 || bc_dist == 0) // distances of zero are invalid
                continue;

            if ((ad_dist == dist) && (bc_dist == dist))
                continue;

            path_a.push_back(id_c);
            path_b.push_back(id_d);

            bfs_ring_check_path(path_a, path_b, all_pairs, offsets, edges);
            return;
        }
    }

    // continuation of the ring could not be found
    path_a.clear();
    path_b.clear();
}

std::vector<std::vector<int>> graph::bfs_rings(const ud_graph_t &graph, int max_ring)
{
    vector<int> offsets = graph.get_int_offsets(),
        edges = graph.get_int_edge_ids();

    // number of vertices
    int nv = static_cast<int>(offsets.size()) - 1;
    if (nv <= 0 || (max_ring != -1 && max_ring < 3))
        return vector<vector<int>>();

    // First get triangles, which are calculated separately
    vector<vector<int> > rings;
    for (int i = 0; i < nv; ++i)
    {
        vector<int> neighbors_a(&edges[offsets[i]], &edges[offsets[i + 1]]);
        for (int j = offsets[i]; j < offsets[i + 1]; ++j)
        {
            int id_b = edges[j];
            if (i >= id_b)
                continue;

            vector<int> neighbors_b(&edges[offsets[id_b]], &edges[offsets[id_b + 1]]),
                intersection;

            set_intersection(
                neighbors_a.begin(), neighbors_a.end(),
                neighbors_b.begin(), neighbors_b.end(),
                back_inserter(intersection)
            );

            for (int id_c : intersection)
                if (id_c > id_b)
                    rings.push_back({i, id_b, id_c});
        }
    }

    if (max_ring == 3)
        return rings;

    // Get all-pairs shortest paths, up to the ring limit. The distance we
    // need to search is only half the size of the maximum ring we are looking for
    int max_dist = (max_ring == -1 ? -1 : max_ring / 2);
    vector<int> all_pairs = bfs_all_pairs(offsets, edges, max_dist);

    for (int i = 0; i < nv; ++i)
    {
        for (int j = i + 1; j < nv; ++j)
        {
            int dist = all_pairs[i * nv + j];

            if (dist >= 2)
            {
                vector<int> path_a(1, i), path_b(1, j);
                bfs_ring_check_path(path_a, path_b, all_pairs, offsets, edges);

                if (!path_a.empty())
                {
                    // since the function only returns one ring for a given pair of vertices,
                    // if two vertices _share_ two or more rings, they get excluded
                    for (size_t k = 0; k < path_a.size(); ++k)
                    {
                        all_pairs[path_a[k] * nv + path_b[k]] = 0;
                        all_pairs[path_b[k] * nv + path_a[k]] = 0;
                    }

                    copy(path_b.begin(), path_b.end(), back_inserter(path_a));
                    rings.push_back(path_a);
                }
            }
        }
    }

    return rings;
}

void graph::write_gplot_graph(std::string filename, const ud_graph_t &graph)
{
    vector<int> offsets = graph.get_int_offsets(),
        edge_ids = graph.get_int_edge_ids();

    int n_vert = static_cast<int>(offsets.size()) - 1;
    if (n_vert <= 1)
        return;

    ofstream outfile(filename.c_str());

    vector<int> adj_map(n_vert * n_vert, 0);
    for (int i = 0; i < n_vert; ++i)
    {
        for (int j = offsets[i]; j < offsets[i + 1]; ++j)
            adj_map[i * n_vert + edge_ids[j]] = std::abs(i - edge_ids[j]);

        for (int j = 0; j < n_vert; ++j)
            outfile << i << '\t' << j << '\t' << adj_map[i * n_vert + j] << '\n';
        outfile << '\n';
    }

    outfile.close();
}
