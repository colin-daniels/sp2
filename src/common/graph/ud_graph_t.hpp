#ifndef SP2_UD_GRAPH_T_HPP
#define SP2_UD_GRAPH_T_HPP

#include "common/util/templates.hpp"
#include <boost/mpi/communicator.hpp>

#include <vector>
#include <cassert>
#include <algorithm>
#include <utility>

namespace sp2 {
namespace graph {

struct ud_edge_t
{
    int id, ///< edge index
        a,  ///< vertex a index
        b;  ///< vertex b index
    operator std::pair<int, int>() const {
        return {a, b};}
};

class ud_graph_t
{
private:
    /// number of vertices
    std::size_t nv;

    /// vector of edge pairs (id_a, id_b). Note that they are
    /// duplicated (a->b and b->a) as well as always sorted
    std::vector<std::pair<int, int>> edge_pairs;

    /// vector to iterators that point to the beginning of the 'set' of
    /// edges for a given vertex, e.g.
    /// (offsets[a])..........(offsets[b])........ =
    /// (   a->b   ), (a->c), (   b->c   ), (b->d)
    std::vector<decltype(edge_pairs)::const_iterator> offsets;

public:
    ud_graph_t();
    explicit ud_graph_t(int num_vertices);
    ud_graph_t(const std::vector<std::pair<int, int>> pairs_in,
        int num_vertices);
    /// note: assumes duplicated edges (a->b and b->a)
    ud_graph_t(const std::vector<int> &offsets_in,
        const std::vector<int> &edge_ids);

    ud_graph_t(ud_graph_t &&other) = default;
    ud_graph_t(const ud_graph_t &other);

    ud_graph_t& operator=(ud_graph_t&& other);
    ud_graph_t& operator=(const ud_graph_t& other);

    bool operator==(const ud_graph_t &other) const;
    bool operator!=(const ud_graph_t &other) const {
        return !(*this == other);}

    struct edge_iterator_t
    {
        typedef int                         difference_type;
        typedef ud_edge_t                   value_type;
        typedef const ud_edge_t * const     pointer;
        typedef const ud_edge_t&            reference;
        typedef std::forward_iterator_tag   iterator_category;

        // std::vector<std::pair<int, int>>::const_iterator
        using pair_iter = decltype(ud_graph_t::offsets)::value_type;

        int index = -1;
        pair_iter iter;

        ud_edge_t operator*() const {
            return {index, iter->first, iter->second};}

        bool operator!=(const edge_iterator_t& o) const {
            return iter != o.iter || index == -1 || o.index == -1;}

        edge_iterator_t operator++() {
            ++iter; ++index; return *this;}
    };

    int add_vertex();
    void add_edge(int a, int b);

    void remove_vertex(int id);
    void remove_edge(int a, int b);

    std::size_t n_vertices() const {return nv;}
    std::size_t n_edges() const {return edge_pairs.size();}

    range_wrap_t<id_iter_t> vertices() const {
        return range_wrap_t<id_iter_t>(0, nv);}

    range_wrap_t<edge_iterator_t> edges() const;
    range_wrap_t<edge_iterator_t> edges(int vert) const;

    int degree(int vert) const;
    bool adjacent(int a, int b) const;

    std::vector<int> get_int_offsets() const;
    std::vector<int> get_int_edge_ids() const;

    std::vector<std::pair<int, int>> get_edge_pairs() const {
        return edge_pairs;}

    std::vector<int> neighbors(int vert) const
    {
        std::vector<int> vert_neighbors;
        for (auto it = offsets[vert]; it != offsets[vert + 1]; ++it)
            vert_neighbors.push_back(it->second);

        return vert_neighbors;
    }

    void bcast(const boost::mpi::communicator &comm, int root);

#ifdef SP2_ENABLE_TESTS
    bool offsets_ok();
    bool edges_are_sorted();
    bool edges_are_duplicated();
#endif // SP2_ENABLE_TESTS

private:
    edge_iterator_t get_iter(std::size_t vertex) const
    {
        return edge_iterator_t{
            static_cast<int>(
                std::distance(edge_pairs.begin(), offsets[vertex])),
            offsets[vertex]
        };
    }

    /// update the pair iterators stored in the offsets vector
    void update_offsets();

    void add_directed(int a, int b);
    void remove_directed(int a, int b);
};

inline range_wrap_t<ud_graph_t::edge_iterator_t>
ud_graph_t::edges() const {
    return {get_iter(0), get_iter(nv)};}

inline range_wrap_t<ud_graph_t::edge_iterator_t>
ud_graph_t::edges(int vert) const {
    return {get_iter(vert), get_iter(vert + 1)};}

inline int ud_graph_t::degree(int vert) const {
    return static_cast<int>(std::distance(offsets[vert], offsets[vert + 1]));}

inline bool ud_graph_t::adjacent(int a, int b) const
{
    return std::binary_search(
        offsets[a], offsets[a + 1],
        std::make_pair(a, b)
    );
}

} // namespace graph
} // namespace sp2

#endif // SP2_UD_GRAPH_T_HPP
