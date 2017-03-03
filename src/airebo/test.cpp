#include <common/math/vec3_t.hpp>
#include <common/graph/ud_graph_t.hpp>
#include "airebo/test.hpp"

#include <cstdint>
#include <templ/algorithm.hpp>
#include <bitset>

struct airebo_atom_t
{
    ptnl::atom_types type;

    sp2::vec3_t pos,
        grad;

    double N_h,
        N_c;


};

class airebo_system_t
{
private:
    std::vector<sp2::vec3_t> position,
        gradient;

    std::vector<ptnl::atom_types> atom_type;

    sp2::graph::ud_graph_t bond_graph;
    std::vector<ptnl::bond_types> bond_type;
    std::vector<double> bond_length;
    std::vector<sp2::vec3_t> bond_vector;


    double total_potential;
public:

    void update();

    template<ptnl::bond_types type>
    void update_bond(sp2::graph::ud_edge_t bond)
    {
        const double r = bond_length[bond.id];

        // the bond cutoff function f^c(r) and its derivative
        auto cutoff = ptnl::cutoff_function<type>(r);

        if (ptnl::val(cutoff) == 0)
            return;

        // the attractive potential V^A(r) and its derivative
        auto attr_ptnl = ptnl::attractive_potential<type>(r);

        // the repulsive potential V^A(r) and its derivative
        auto rep_ptnl = ptnl::repulsive_potential<type>(r);

        // sqrt(x^2 + y^2 + z^2)
        // d/dx = x / sqrt(x^2 + y^2 + z^2)

        // add the repulsive potential to the total potential now,
        // since no other terms are needed
        total_potential += ptnl::val(cutoff) * ptnl::val(rep_ptnl);

//        gradient[bond.a] += cutoff.deriv() * rep_ptnl.value() +
//                            cutoff.value() * rep_ptnl.deriv();

//        bond_dir_force[j] += cutoff[j * 2 + 1] * rep_p + cutoff[j * 2] * d_rep_p;
    }
};

void airebo_system_t::update()
{
    bond_type.clear();

    for (auto bond : bond_graph.edges())
    {
        const int id_a = bond.a,
            id_b = bond.b;

        ptnl::bond_types type = ptnl::get_bond_type(
            atom_type[id_a], atom_type[id_b]);

        switch (type)
        {
        case ptnl::bond_types::CC:
            update_bond<ptnl::bond_types::CC>(bond);
            break;
        case ptnl::bond_types::HC:
            update_bond<ptnl::bond_types::HC>(bond);
            break;
        case ptnl::bond_types::HH:
            update_bond<ptnl::bond_types::HH>(bond);
            break;
        }
    }
}

constexpr std::uint32_t elem_hash(const char* name)
{
    return  static_cast<std::uint32_t>(name[0]) |
        static_cast<std::uint32_t>(name[1]) << CHAR_BIT;
}

enum class test_enum : std::uint32_t
{
    C = elem_hash("C")
};

void testfn()
{
    double val, deriv;

    std::tie(
        val,
        deriv
    ) = ptnl::attractive_potential<ptnl::bond_types::HH>(10);

    const char *input = "He";

    switch(elem_hash(input))
    {
    case elem_hash("H"):
        return;
    }
}