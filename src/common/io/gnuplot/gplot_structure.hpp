#ifndef SP2_GPLOT_STRUCTURE_HPP
#define SP2_GPLOT_STRUCTURE_HPP

#include <common/structure_t.hpp>
#include <common/graph/ud_graph_t.hpp>
#include <string>
#include <functional>
#include <iostream>

namespace sp2 {
namespace io {

void draw_top_down(std::ostream &output,
    structure_t structure,
    std::pair<double, double> x_bounds,
    std::pair<double, double> y_bounds,
    const graph::ud_graph_t &bond_graph,
    double atom_radius = 0.25, double bond_radius = 0.2,
    std::function<void(int, std::ostream&)> foreach_atom
        = [](int, std::ostream&){}
);

} // namespace io
} // namespace sp2

#endif // SP2_GPLOT_STRUCTURE_HPP
