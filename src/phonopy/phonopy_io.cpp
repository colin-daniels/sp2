#include "phonopy_io.hpp"
#include "yaml-cpp/yaml.h"

#include <map>
#include <iostream>
#include <fstream>
#include <common/io/gnuplot/gplot_structure.hpp>
#include <common/math/vec3_util.hpp>
#include <airebo/system_control_t.hpp>
#include <common/util/modeling.hpp>
#include <common/atom_types.hpp>
#include <common/math/mat3x3_t.hpp>

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

inline sp2::vec3_t node_to_v3(const YAML::Node &node)
{
    return {node[0].as<double>(), node[1].as<double>(), node[2].as<double>()};
}

sp2::structure_t read_phonopy_structure(YAML::Node &root)
{
    sp2::structure_t structure;

    int n_atoms = root["natom"].as<int>();

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            structure.lattice[i][j] = root["lattice"][i][j].as<double>();

    sp2::mat3x3_t lat(structure.lattice);
    // transpose because lattice vectors are rows, need columns to convert to
    // cartesian coordinates
    lat.transpose();

    for (const auto &atom : root["points"])
    {
        auto symbol = atom["symbol"].as<std::string>();
        structure.types.push_back(
            symbol == "H" ? sp2::atom_type::HYDROGEN :
            symbol == "C" ? sp2::atom_type::CARBON :
                throw std::runtime_error("Bad atom symbol " +
                    symbol + "in phonopy file.")
        );

        // transform to cartesian and store coordinates
        for (auto coord : lat * node_to_v3(atom["coordinates"]))
            structure.positions.push_back(coord);
    }

    if (n_atoms != structure.types.size() ||
        n_atoms * 3 != structure.positions.size())
        throw std::runtime_error("Couldn't read phonopy structure.");

    return structure;
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
    std::pair<double, std::vector<sp2::vec3_t>> mode)
{
    constexpr double atom_radius = 0.2,
        bond_radius = atom_radius,
        arrow_radius = 2 * atom_radius,
        amplification = 3 * arrow_radius;

    // need to center the structure first to avoid cutting it in two with
    // the lattice
    structure = util::center_by_avg(structure);

    // remove hydrogen atoms
    auto pos = dtov3(structure.positions);
    for (int i = 0; i < structure.types.size(); ++i)
    {
        if (structure.types[i] != atom_type::HYDROGEN)
            continue;

        pos.erase(pos.begin() + i);
        structure.types.erase(structure.types.begin() + i);
        mode.second.erase(mode.second.begin() + i);

        i -= 1;
    }
    structure.positions = sp2::v3tod(pos);

    // rescale eigenvectors by largest magnitude
    double mag_max = 0;
    for (auto &eig : mode.second)
        mag_max = std::max(mag_max, eig.mag_sq());

    mag_max = std::sqrt(mag_max);
    for (auto &eig : mode.second)
        eig /= mag_max;

    std::ofstream outfile(filename);

    // set up the plot area
    auto bounds = get_bounds({
        vec3_t(structure.lattice[0]),
        vec3_t(structure.lattice[1]),
        vec3_t(structure.lattice[2])
    });

    outfile << "set size ratio -1;\n"
            << "set xrange ["
                << bounds.first.x << ":" << bounds.second.x
            << "];\n"
            << "set yrange ["
                << bounds.first.y << ":" << bounds.second.y
            << "];\n"
            << "set style arrow 1 head filled size "
                    << arrow_radius
            << ", 30, 0 fixed lw 3 lc rgb \"#1f78b4\";\n";

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
            auto eigenvector = mode.second[atom_id] * amplification
                    + mode.second[atom_id].unit_vector() * arrow_radius,
                arrow_start = pos[atom_id]
                    + eigenvector.unit_vector() * atom_radius;

            if (eigenvector.mag() < 1.5 * arrow_radius)
                return;

            // draw an arrow representing the mode eigenvector
            oss << "set arrow from "
                    << arrow_start.x << ", " << arrow_start.y
                << " rto "
                    << eigenvector.x << ", " << eigenvector.y
                << " as 1;\n";
    });
}

void sp2::phonopy::write_force_sets(
    std::string filename,
    const std::pair<structure_t, std::vector<std::pair<int, vec3_t>>>
        &displacements,
    const std::vector<std::vector<vec3_t>> &forces
)
{
    // From the Phonopy Manual
    //
    // This file gives sets of forces in supercells with finite atomic
    // displacements. Each supercell involves one displaced atom. The first
    // line is the number of atoms in supercell. The second line gives number
    // of calculated supercells with displacements. Below the lines, sets of
    // forces with displacements are written. In each set, firstly the atom
    // number in supercell is written. Secondary, the atomic displacement in
    // Cartesian coordinates is written. Below the displacement line, atomic
    // forces in Cartesian coordinates are successively written. This is
    // repeated for the set of displacements. Blank likes are simply ignored.

    if (displacements.second.size() != forces.size())
        throw std::runtime_error("Number of displacements and forces "
            "calculated do not match.");

    std::ofstream outfile(filename);
    outfile.precision(std::numeric_limits<double>::digits10);

    // number of atoms
    outfile << displacements.first.types.size() << '\n'
            // number of displacements
            << displacements.second.size() << '\n'
            << '\n';

    for (std::size_t i = 0; i < forces.size(); ++i)
    {
        const auto &disp = displacements.second[i];

        // atom id, note: phonopy 1-indexes their atoms
        outfile << disp.first + 1 << '\n'
                // cartesian displacement for the given atom
                << disp.second.x << ' '
                    << disp.second.y << ' '
                    << disp.second.z << '\n';

        // forces on all atoms
        for (auto &v : forces[i])
            outfile << v.x << ' ' << v.y << ' ' << v.z << '\n';

        // blank line for easy reading
        outfile << '\n';
    }
}

std::pair<sp2::structure_t, std::vector<std::pair<int, sp2::vec3_t>>>
    sp2::phonopy::read_displacements(std::string filename)
{
    YAML::Node node = YAML::LoadFile(filename);

    structure_t structure = read_phonopy_structure(node);
    std::vector<std::pair<int, vec3_t>> displacements;

    for (const auto &disp : node["displacements"])
    {
        displacements.emplace_back(
            // atom id (phonopy 1-indexes their atoms)
            disp["atom"].as<int>() - 1,
            // convert cartesian displacement array into a vector
            node_to_v3(disp["displacement"])
        );
    }

    return std::make_pair(structure, displacements);
}
