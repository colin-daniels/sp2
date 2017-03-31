#include "phonopy_io.hpp"
#include "yaml-cpp/yaml.h"

#include <map>
#include <iostream>
#include <fstream>
#include <common/io/gnuplot/gplot_structure.hpp>
#include <common/math/vec3_util.hpp>
#include <airebo/system_control_t.hpp>
#include <common/util/modeling.hpp>

void print_type(const YAML::Node &node)
{
    std::string type = "";
    switch (node.Type()) {
    case YAML::NodeType::Null:      type = "Null"; break;
    case YAML::NodeType::Scalar:    type = "Scalar"; break;
    case YAML::NodeType::Sequence:  type = "Sequence"; break;
    case YAML::NodeType::Map:       type = "Map"; break;
    case YAML::NodeType::Undefined: type = "Undefined"; break;
    }

    std::cerr << "Node type: " << type << std::endl;
}

std::vector<std::string> sp2::phonopy::read_irreps(std::string filename)
{
    YAML::Node node = YAML::LoadFile(filename);

    if (!node || !node["normal_modes"])
        return {};


    std::map<int, std::string> irrep_labels;
    for (const auto &elem : node["normal_modes"])
    {
        if (!elem["band_indices"])
            return {};

        std::string label = elem["ir_label"] ?
            elem["ir_label"].as<std::string>() : "Unknown";

        for (const auto &idx_node : elem["band_indices"])
        {
            int idx = idx_node.as<int>();
            irrep_labels[idx] = label;
        }
    }

    std::vector<std::string> output;
    for (auto &elem : irrep_labels)
        output.push_back(elem.second);

    return output;
}

void sp2::phonopy::draw_normal_mode(std::string filename,
    sp2::structure_t structure,
    const std::pair<double, std::vector<sp2::vec3_t>> &mode)
{
    constexpr double atom_radius = 0.25,
        bond_radius = 0.2,
        amplification = 5.0;

    // need to center the structure first to avoid cutting it in two with
    // the lattice
    structure = util::center_by_avg(structure);

    auto pos = dtov3(structure.positions);

    std::ofstream outfile(filename);

    // set up the plot area
    auto bounds = get_bounds({
        vec3_t(structure.lattice[0]),
        vec3_t(structure.lattice[1]),
        vec3_t(structure.lattice[2])
    });

    outfile << "set xrange ["
                << bounds.first.x << ":" << bounds.second.x
            << "];\n"
            << "set yrange ["
                << bounds.first.y << ":" << bounds.second.y
            << "];\n";

    // set structure lattice to zero so bonds don't wrap
    for (int i = 0; i < 9; ++i)
        structure.lattice[i / 3][i % 3] = 0;

    io::draw_top_down(outfile, structure,
        {bounds.first.x, bounds.second.x},
        {bounds.first.y, bounds.second.y},
        airebo::system_control_t(structure).get_bond_control().get_graph(),
        atom_radius, bond_radius,
        // for each atom
        [&](int atom_id, std::ostream &oss) {
            auto eigenvector = mode.second[atom_id] * amplification,
                arrow_start = pos[atom_id]
                    + eigenvector.unit_vector() * atom_radius * 2;

            // draw an arrow representing the mode eigenvector
            oss << "set arrow from "
                    << arrow_start.x << ", " << arrow_start.y
                << " rto "
                    << eigenvector.x << ", " << eigenvector.y
                << ";\n";
    });
}
