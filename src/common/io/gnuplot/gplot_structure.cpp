#include "gplot_structure.hpp"
#include <sstream>
#include <common/math/vec3_util.hpp>

inline void init_object_num(std::ostream &output)
{
    // initialize object number if not present
    output << "# Initialize an object number\n"
           << "if (!exists(\"object_num\")) object_num = 1;\n";
}

inline void output_object_num(std::ostream &output, int num)
{
    output << "(object_num + " << num << ")";
}

inline void update_object_num(std::ostream &output, int delta)
{
    output << "object_num = object_num + " << delta << ";\n";
}

void sp2::io::draw_top_down(std::ostream &output,
    sp2::structure_t structure,
    std::pair<double, double> x_bounds,
    std::pair<double, double> y_bounds,
    const sp2::graph::ud_graph_t &bond_graph,
    double atom_radius, double bond_radius,
    std::function<void(int, std::ostream&)> foreach_atom)
{
    init_object_num(output);
    int current_object_num = 0;

    auto positions = sp2::dtov3(structure.positions);
    auto &types = structure.types;

    auto in_bounds = [&](const vec3_t &pos) {
        return x_bounds.first < pos.x && pos.x < x_bounds.second &&
            y_bounds.first < pos.y && pos.y < y_bounds.second;
    };

    auto clamp_pos = [&](vec3_t &pos) {
        pos.x = std::max(x_bounds.first, std::min(x_bounds.second, pos.x));
        pos.y = std::max(y_bounds.first, std::min(y_bounds.second, pos.y));
    };

    // draw bonds
    for (const auto &edge : bond_graph.edges())
    {
        vec3_t atom_a = positions[edge.a],
            atom_b = positions[edge.b];

        if (!in_bounds(atom_a) && !in_bounds(atom_b))
            continue;

        vec3_t normal = bond_radius * unit_normal_to(
            atom_b - atom_a, vec3_t(0, 0, 1));

        output << "set object ";
        output_object_num(output, current_object_num);
        output << " polygon from ";

        bool first = true;
        for (auto pos : {
            atom_a + normal, atom_b + normal,
            atom_b - normal, atom_a - normal
        })
        {
            if (!first)
                output << " to ";
            else
                first = false;

            // make sure we don't draw outside the lines
            clamp_pos(pos);

            output << pos.x << ", " << pos.y;
        }

        output << '\n';
        current_object_num++;
    }

    for (int i = 0; i < positions.size(); ++i)
    {
        auto &pos = positions[i];
        if (!in_bounds(pos))
            continue;

        // draw a circle at the atom
        output << "set object ";
        output_object_num(output, current_object_num);
        output << " circle at " << pos.x << ", " << pos.y
               << " size " << atom_radius << "\n";

        foreach_atom(i, output);
        current_object_num++;
    }

    update_object_num(output, current_object_num);
}
