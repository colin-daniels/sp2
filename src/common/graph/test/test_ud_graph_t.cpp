#include "common/graph/ud_graph_t.hpp"

#include <gtest/gtest.h>

using namespace std;
using namespace sp2::graph;

struct simple_graph_t
{
    vector<int> offsets,
        edge_ids;
};

vector<pair<int, int>> gen_graph_pairs(int nv = 1000, int ne = -1)
{
    if (ne < 0)
        ne = nv * 8;

    vector<pair<int, int>> pairs;
    for (int i = 0; i < ne; ++i)
    {
        auto p = make_pair(rand() % nv, rand() % nv);
        pairs.push_back(p);
        pairs.push_back({p.second, p.first});
    }

    // deduplicate
    sort(pairs.begin(), pairs.end());
    pairs.erase(unique(pairs.begin(), pairs.end()), pairs.end());

    return pairs;
}

simple_graph_t gen_graph_simple(int nv = 1000, int ne = -1)
{
    if (ne < 0)
        ne = nv * 8;

    auto pairs = gen_graph_pairs(nv, ne);

    vector<int> edge_ids,
        offsets(nv, 0);
    for (auto p : pairs)
    {
        edge_ids.push_back(p.second);
        offsets[p.first]++;
    }
    offsets.insert(offsets.begin(), 0);

    for (size_t i = 1; i < offsets.size(); ++i)
        offsets[i] += offsets[i - 1];

    return  {offsets, edge_ids};
}

void validate_graph(ud_graph_t &graph)
{
    ASSERT_TRUE(graph.offsets_ok());
    ASSERT_TRUE(graph.edges_are_sorted());
    ASSERT_TRUE(graph.edges_are_duplicated());
}

TEST(test_graph, edgeless)
{
    ud_graph_t graph;

    EXPECT_EQ(0, graph.n_vertices());
    EXPECT_EQ(0, graph.n_edges());
    validate_graph(graph);

    for (int i = 0; i < 100; ++i)
        graph.add_vertex();

    EXPECT_EQ(100, graph.n_vertices());
    EXPECT_EQ(ud_graph_t(100), graph);
    validate_graph(graph);
}

TEST(test_graph, offset_construct)
{
    int n_vert = 1000;
    auto ref = gen_graph_simple(n_vert);

    // test constructor
    ud_graph_t graph(ref.offsets, ref.edge_ids);

    for (size_t i = 0; i + 1 < ref.offsets.size(); ++i)
    {
        ASSERT_EQ(ref.offsets[i + 1] - ref.offsets[i], graph.degree(i));
        for (int j = ref.offsets[i]; j < ref.offsets[i + 1]; ++j)
        {
            ASSERT_TRUE(graph.adjacent(i, ref.edge_ids[j]));
            ASSERT_TRUE(graph.adjacent(ref.edge_ids[j], i));
        }
    }

    validate_graph(graph);
    EXPECT_EQ(ref.offsets, graph.get_int_offsets());
    EXPECT_EQ(ref.edge_ids, graph.get_int_edge_ids());

    EXPECT_EQ(n_vert, graph.n_vertices());
    EXPECT_EQ((int)ref.edge_ids.size(), graph.n_edges());
}

TEST(test_graph, equality_assignment)
{
    auto ref_1 = gen_graph_simple(),
        ref_2 = gen_graph_simple();

    ud_graph_t graph_a(ref_1.offsets, ref_1.edge_ids),
        graph_b(ref_2.offsets, ref_2.edge_ids);

    EXPECT_NE(graph_a, graph_b);

    EXPECT_EQ(graph_a, graph_a);
    EXPECT_EQ(ud_graph_t(),   ud_graph_t());
    EXPECT_NE(ud_graph_t(10), ud_graph_t());

    ud_graph_t graph_c(graph_a);
    EXPECT_EQ(graph_a, graph_c);

    graph_c = graph_b;
    EXPECT_EQ(graph_b, graph_c);

    graph_a = ud_graph_t(graph_c);
    EXPECT_EQ(graph_a, graph_b);

    validate_graph(graph_a);
    validate_graph(graph_b);
    validate_graph(graph_c);
}

TEST(test_graph, pair_construct)
{
    int nv = 1000;
    auto pairs = gen_graph_pairs(nv);

    // duplicated/sorted pairs
    ud_graph_t graph_a(pairs, nv);
    validate_graph(graph_a);
    EXPECT_EQ(pairs.size(), graph_a.n_edges());

    // duplicated/randomized pairs
    random_shuffle(pairs.begin(), pairs.end());

    ud_graph_t graph_b(pairs, nv);
    validate_graph(graph_b);
    EXPECT_EQ(pairs.size(), graph_b.n_edges());

    // unique/randomized pairs
    vector<pair<int, int>> pairs_single;
    copy_if(pairs.begin(), pairs.end(), back_inserter(pairs_single),
        [](pair<int, int> p) {return p.first >= p.second;});

    ud_graph_t graph_c(pairs_single, nv);
    validate_graph(graph_c);
    EXPECT_EQ(pairs.size(), graph_c.n_edges());

    EXPECT_EQ(graph_a, graph_b);
    EXPECT_EQ(graph_b, graph_c);
    EXPECT_EQ(graph_c, graph_a);
}

TEST(test_graph, edge_iterators)
{
    int nv = 1000;
    auto pairs = gen_graph_pairs(nv);

    ud_graph_t graph(pairs, nv);

    // test entire graph iterator
    int idx_outer = 0;
    for (ud_edge_t edge : graph.edges())
    {
        ASSERT_EQ(idx_outer, edge.id);
        ASSERT_EQ(pairs[idx_outer].first, edge.a);
        ASSERT_EQ(pairs[idx_outer].second, edge.b);

        ++idx_outer;
    }

    ASSERT_EQ(idx_outer, (int)pairs.size());

    vector<int> offsets = graph.get_int_offsets(),
        edge_ids = graph.get_int_edge_ids();

    // test iterator for individual vertices
    for (size_t i = 0; i + 1 < offsets.size(); ++i)
    {
        int idx_inner = offsets[i];
        for (ud_edge_t edge : graph.edges(i))
        {
            ASSERT_EQ(idx_inner, edge.id);
            ASSERT_EQ(i, edge.a);
            ASSERT_EQ(edge_ids[idx_inner], edge.b);

            idx_inner += 1;
        }
        ASSERT_EQ(idx_inner, offsets[i + 1]);
    }
}

vector<pair<int, int>> a_not_in_b(ud_graph_t &graph_a, ud_graph_t &graph_b)
{
    vector<pair<int, int>> edges,
        edges_a, edges_b;

    for (auto p : graph_a.get_edge_pairs())
        if (p.first <= p.second) edges_a.push_back(p);

    for (auto p : graph_b.get_edge_pairs())
        if (p.first <= p.second) edges_b.push_back(p);

    set_difference(edges_a.begin(), edges_a.end(),
        edges_b.begin(), edges_b.end(),
        back_inserter(edges));

    return edges;
}

// TODO: Fix add/rem_transform

//TEST(test_graph, rem_transform)
//{
//    int nv = 100, ne = nv * nv; // dense
//    ud_graph_t graph_a(gen_graph_pairs(nv, ne), nv),
//        graph_b(gen_graph_pairs(nv, ne), nv);
//
//    vector<pair<int, int>> to_remove_a = a_not_in_b(graph_a, graph_b),
//        to_remove_b = a_not_in_b(graph_b, graph_a);
//
//    for (auto p : to_remove_a)
//        graph_a.remove_edge(p.first, p.second);
//
//    for (auto p : to_remove_b)
//        graph_b.remove_edge(p.first, p.second);
//
//    validate_graph(graph_a);
//    validate_graph(graph_b);
//
//    ASSERT_EQ(graph_a, graph_b);
//}
//
//TEST(test_graph, add_transform)
//{
//    int nv = 100, ne = nv * nv; // dense
//    ud_graph_t graph_a(gen_graph_pairs(nv, ne), nv),
//        graph_b(gen_graph_pairs(nv, ne), nv);
//
//    vector<pair<int, int>> to_add_a = a_not_in_b(graph_b, graph_a),
//        to_add_b = a_not_in_b(graph_a, graph_b);
//
//    for (auto p : to_add_a)
//        graph_a.add_edge(p.first, p.second);
//
//    for (auto p : to_add_b)
//        graph_b.add_edge(p.first, p.second);
//
//    validate_graph(graph_a);
//    validate_graph(graph_b);
//
//    ASSERT_EQ(graph_a, graph_b);
//}
//
//TEST(test_graph, add_rem_transform)
//{
//    int nv = 100, ne = nv * nv; // dense
//    ud_graph_t graph_a(gen_graph_pairs(nv, ne), nv),
//        graph_b(gen_graph_pairs(nv, ne), nv);
//
//    vector<pair<int, int>> to_remove = a_not_in_b(graph_a, graph_b),
//        to_add = a_not_in_b(graph_b, graph_a);
//
//    for (auto p : to_add)
//        graph_a.add_edge(p.first, p.second);
//
//    for (auto p : to_remove)
//        graph_a.remove_edge(p.first, p.second);
//
//    validate_graph(graph_a);
//
//    ASSERT_EQ(graph_a, graph_b);
//}