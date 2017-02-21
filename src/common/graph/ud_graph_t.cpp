#include "common/graph/ud_graph_t.hpp"

#include <boost/mpi.hpp>
#include <boost/serialization/utility.hpp>
#include <cassert>

using namespace std;
using namespace sp2;
namespace mpi = boost::mpi;

vector<pair<int, int>> make_pairlist(const vector<int> &offsets,
    const vector<int> &edge_ids)
{
    if (offsets.back() != static_cast<int>(edge_ids.size()))
        throw invalid_argument("invalid argument passed to make_pairlist, "
            "offsets last element does not equal number of edges.");

    vector<pair<int, int>> pairs;
    for (size_t i = 0; i + 1 < offsets.size(); ++i)
        for (int j = offsets[i]; j < offsets[i + 1]; ++j)
            pairs.emplace_back(i, edge_ids[j]);

    auto bad = any_of(pairs.begin(), pairs.end(), [](auto p) {
        if (p.first < 0 || p.second < 0)
        {
            cout << "bad pair: " << p.first << ' ' << p.second << endl;
            return true;
        }
        return false;
    });

    if (bad)
    {
        cout << "bad pairlist" << endl;

        bad = any_of(edge_ids.begin(), edge_ids.end(), [&](auto i) {
            if (i < 0 || i >= offsets.size())
            {
                cout << i << endl;
                return true;
            }
            return false;
        });

        if (bad)
            cout << "bad edge ids" << endl;

        bad = any_of(offsets.begin(), offsets.end(), [&](auto i) {
            return i < 0 || i > edge_ids.size();
        });

        if (bad)
            cout << "bad offsets" << endl;

        abort();
    }

    return pairs;
}

graph::ud_graph_t::ud_graph_t() :
    nv(0), offsets(1, edge_pairs.begin()) {}

graph::ud_graph_t::ud_graph_t(int num_vertices) :
    nv(num_vertices),
    offsets(nv + 1, edge_pairs.begin()) {}

graph::ud_graph_t::ud_graph_t(const std::vector<std::pair<int, int>> pairs_in,
    int num_vertices) : nv(num_vertices)
{
    if (pairs_in.empty())
    {
        offsets.resize(nv + 1, edge_pairs.begin());
        return;
    }

    // duplicate the pairs
    edge_pairs = pairs_in;
    for (auto &p : pairs_in)
        edge_pairs.push_back({p.second, p.first});

    sort(edge_pairs.begin(), edge_pairs.end());

    // deduplicate
    auto last = unique(edge_pairs.begin(), edge_pairs.end());
    edge_pairs.erase(last, edge_pairs.end());

    // populate offsets
    update_offsets();
}

graph::ud_graph_t::ud_graph_t(const std::vector<int> &offsets_in,
    const std::vector<int> &edge_ids) :
    ud_graph_t(make_pairlist(offsets_in, edge_ids),
        offsets_in.size() - 1) {}



//graph::ud_graph_t::ud_graph_t(ud_graph_t &&other) {*this = other;}
graph::ud_graph_t::ud_graph_t(const ud_graph_t &other) :
    nv(other.nv), edge_pairs(other.edge_pairs) {
    update_offsets();}

graph::ud_graph_t& graph::ud_graph_t::operator=(ud_graph_t&& other)
{
    nv = other.nv;
    edge_pairs = move(other.edge_pairs);
    offsets = move(other.offsets);
    return *this;
}

graph::ud_graph_t& graph::ud_graph_t::operator=(const ud_graph_t& other) {
    return *this = ud_graph_t(other);}

bool graph::ud_graph_t::operator==(const ud_graph_t &other) const {
    return (nv == other.nv && edge_pairs == other.edge_pairs);}

void graph::ud_graph_t::bcast(const boost::mpi::communicator &comm, int root)
{
    // send vertices
    mpi::broadcast(comm, nv, root);
    mpi::broadcast(comm, edge_pairs, root);

    // regenerate offsets
    update_offsets();
}

std::vector<int> graph::ud_graph_t::get_int_offsets() const
{
    std::vector<int> result;
    for (auto &it : offsets)
        result.push_back(std::distance(edge_pairs.begin(), it));

    return result;
}

std::vector<int> graph::ud_graph_t::get_int_edge_ids() const
{
    std::vector<int> result;
    for (auto p : edge_pairs)
        result.push_back(p.second);
    return result;
}

void graph::ud_graph_t::update_offsets()
{
    offsets.clear();
    for (auto it = edge_pairs.begin(); it != edge_pairs.end(); ++it)
        if (static_cast<int>(offsets.size()) < it->first + 1)
            offsets.resize(it->first + 1, it);

    offsets.resize(nv + 1, edge_pairs.end());
}

int graph::ud_graph_t::add_vertex()
{
    offsets.push_back(edge_pairs.end());
    return nv++;
}

void graph::ud_graph_t::remove_vertex(int id)
{
    auto first = edge_pairs.begin();
    for (auto i = first + 1; i != edge_pairs.end(); ++i)
    {
        // shift down edge indices
        if (i->second > id) i->second--;
        if (i->first  > id) i->first--;

        // effectively overwrite edges connected to id by shifting
        // values later in the vector down one
        if (i->second != id && i->first != id)
            *first++ = std::move(*i);
    }

    nv -= 1;
    edge_pairs.erase(first, edge_pairs.end());
    update_offsets();
}

void graph::ud_graph_t::add_edge(int a, int b)
{
    add_directed(a, b);
    add_directed(b, a);
}

void graph::ud_graph_t::remove_edge(int a, int b)
{
    remove_directed(a, b);
    remove_directed(b, a);
}

void graph::ud_graph_t::add_directed(int a, int b)
{
    // whether or not the new edge is going to invalidate our iterators
    bool iterators_invalidated = edge_pairs.size() == edge_pairs.capacity();

    // insert it so that edges are still sorted
    auto new_edge = std::make_pair(a, b);
    edge_pairs.insert(
        std::upper_bound(offsets[a], offsets[a + 1], new_edge), new_edge);

    // update offsets differently depending on if we had to resize
    if (iterators_invalidated)
        update_offsets();
    else
    {
        for (size_t i = a + 1; i < offsets.size(); ++i)
            ++offsets[i];
    }
}

void graph::ud_graph_t::remove_directed(int a, int b)
{
    // find the given edge from a -> b
    auto pos = std::lower_bound(
        offsets[a], offsets[a + 1],
        std::make_pair(a, b)
    );

    // failed to find edge
    if (pos == offsets[a + 1])
        return;

    edge_pairs.erase(pos);

    // since removal wont invalidate iterators for vector
    // we just update all offsets past the vertex (a)
    for (std::size_t i = a + 1; i < offsets.size(); ++i)
        --offsets[i];
}

#ifdef SP2_ENABLE_TESTS
bool graph::ud_graph_t::offsets_ok()
{
    for (int i = 0; i < nv; ++i)
    {
        for (auto it = offsets[i]; it != offsets[i + 1]; ++it)
            if (it->first != i)
                return false;

        if (offsets[i] != edge_pairs.begin() &&
            (offsets[i] - 1)->first == i)
            return false;

        if (offsets[i + 1] != edge_pairs.end() &&
            offsets[i + 1]->first == i)
            return false;

    }
    return true;
}

bool graph::ud_graph_t::edges_are_sorted()
{
    return is_sorted(edge_pairs.begin(), edge_pairs.end());
}

bool graph::ud_graph_t::edges_are_duplicated()
{
    return all_of(edge_pairs.begin(), edge_pairs.end(),
        [&](pair<int, int> p) {
            return adjacent(p.second, p.first) &&
                   adjacent(p.first, p.second);
    });
}
#endif // SP2_ENABLE_TESTS

