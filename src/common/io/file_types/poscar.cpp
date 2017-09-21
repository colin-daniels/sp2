//
// Created by cc on 9/20/16.
//

#include "common/io/file_types/poscar.hpp"
#include "common/io/util.hpp"
#include "common/neighbor/utility_functions.hpp"

#include <string>
#include <fstream>
#include <common/math/vec3_t.hpp>

#include <variant>
#include <optional>
#include <src/common/atom_types.hpp>

using namespace sp2;

std::variant<sp2::structure_t, std::string>
    read_poscar_impl(const std::string &filename)
{
    // to be read into from file
    sp2::structure_t input;

    std::ifstream infile(filename);
    if (!infile) { return "could not open"; }

    // read next line from the file, return a null option if it fails
    auto next_line = [&infile]() -> std::optional<std::string> {
        std::string line;
        if (getline(infile, line))
            return line;

        return std::nullopt;
    };

    // first line is a comment
    if (!next_line()) { return "could not read first line"; }

////////////////////////////////////////////////////////////////////////////////
// Lattice [lines 2-5]                                                        //
////////////////////////////////////////////////////////////////////////////////

    // second line is the lattice scaling factor
    double lattice_scale;
    try {
        lattice_scale = std::stod(
            next_line().value()
        );
    } catch (std::exception&) {
        return "could not read/parse lattice scaling factor";
    }

    // third, fourth, and fifth lines are the lattice vectors as row vectors
    for (auto &vec : input.lattice)
    {
        auto line = next_line();
        if (!line) { return "could not read lattice vectors"; }

        auto split = io::split_line(*line);
        if (split.size() != 3) { return "could not parse lattice vectors"; }

        try {
            for (int j = 0; j < 3; ++j)
                vec[j] = lattice_scale * std::stod(split[j]);
        } catch (std::exception&) {
            return "could not parse lattice vectors";
        }
    }

////////////////////////////////////////////////////////////////////////////////
// Atom types/number of atoms                                                 //
////////////////////////////////////////////////////////////////////////////////

    // sixth line can either be the atom numbers or the atom types
    {
        auto line = next_line();
        if (!line) { return "could not read atom types/numbers"; }

        // see if the next line is the atom numbers, or atom types
        bool has_atom_types = std::all_of(
            line->begin(), line->end(), [](char c) {
            return std::isspace(c) || std::isalpha(c);
        });

        std::vector<atom_types> type_list;
        if (has_atom_types)
        {
            // read atom types from strings, note: case insensitive
            for (const std::string &type_str : io::split_line(*line))
                type_list.push_back(enum_from_str<atom_types, true>(type_str));

            // if this line was atom types, we need to read the next line to
            // get the atom numbers
            line = next_line();
            if (!line)
                return "could not read atom numbers";
        }

        // read atom numbers
        std::vector<std::size_t> atom_type_cutoffs = {0};
        try {
            io::trim(*line, " \t");
            for (const auto &num_str : io::split_line(*line))
            {
                atom_type_cutoffs.push_back(
                    atom_type_cutoffs.back() + std::stoi(num_str)
                );
            }
        } catch (std::exception&) {
            return "could not parse atom numbers";
        }

        // atom numbers are required
        std::size_t n_types = atom_type_cutoffs.size() - 1;
        if (n_types == 0)
            return "missing number of atoms in file";

        // default for types if none specified
        if (!has_atom_types)
        {
            if (n_types == 1)
                type_list = {atom_types::C};
            else if (n_types == 2)
                type_list = {atom_types::H, atom_types::C};
            else
                return "no atom types specified, defaults not acceptable";

            std::cout << "Warning: used defaults for non-specified atom types "
                      "in POSCAR file \'" << filename << "\'.\n";
        }
        else if (n_types != type_list.size())
            return "mismatch in number of atom types and list of atom numbers";

        // put together actual atom type vector
        for (std::size_t i = 0; i < type_list.size(); ++i)
            input.types.resize(atom_type_cutoffs[i + 1], type_list[i]);

    } // scope

    // total number of atoms is the same as the final atom type 'cutoff'
    auto total_atoms = input.types.size();

////////////////////////////////////////////////////////////////////////////////
// Coordinates and selective dynamics tags                                    //
////////////////////////////////////////////////////////////////////////////////

    bool selective_dyn = false,
         direct_coords = false;

    // the next one or two lines specify selective dynamics and the coordinate
    // specification
    {
        auto line = next_line();
        if (!line || line->empty())
            return "could not read coordinate specifier (direct/cartesian) or "
                "selective dynamics specifier";

        // this line defines selective dynamics if it starts with S
        if (std::toupper(line->front()) == 'S')
        {
            selective_dyn = true;

            // since this line was the selective dynamics specifier, we need
            // to read in the next one to get the coordinate specifier
            line = next_line();
            if (!line || line->empty())
                return "could not read coordinate specifier (direct/cartesian)";
        }

        direct_coords = std::toupper(line->front()) != 'C' &&
                        std::toupper(line->front()) != 'K';
    } // scope

////////////////////////////////////////////////////////////////////////////////
// Read in atom coordinates                                                   //
////////////////////////////////////////////////////////////////////////////////

    auto lattice_transform = mat3x3_t{input.lattice}.transposed();
    for (std::size_t i = 0; i < total_atoms; ++i)
    {
        auto line = next_line();
        if (!line)
            return "failed to read coordinates for atom number "
                   + std::to_string(i);

        // error message
        const auto parse_fail = "failed to parse coordinates for atom number "
               + std::to_string(i) + ", text: \"" + *line + "\"";

        auto split = io::split_line(*line);
        if (split.size() != 3 || (split.size() != 6 && selective_dyn))
            return parse_fail;

        try {
            vec3_t pos = {
                std::stod(split[0]),
                std::stod(split[1]),
                std::stod(split[2])
            };

            // for direct coords, x = L^T x
            if (direct_coords)
                pos = lattice_transform * pos;

            input.positions.push_back(pos);

            // TODO: Selective dynamics

        } catch (std::exception&) {
            return parse_fail;
        }
    }

    return input;
}


bool io::read_poscar(std::string filename, structure_t &output)
{
    auto result = read_poscar_impl(filename);

    if (std::holds_alternative<std::string>(result))
    {
        // an error occured
        std::cerr << "Error when reading POSCAR file \""
                  << filename << "\": " << std::get<std::string>(result)
                  << std::endl;

        return false;
    }
    else if (std::holds_alternative<structure_t>(result))
    {
        // seems O.K.
        output = std::get<structure_t>(result);
        return true;
    }
    else
    {
        // variant isn't holding anything
        std::cerr << "Variant issue in read_poscar.\n";
        return false;
    }
}

bool io::write_poscar(std::string filename, const structure_t &input)
{
    if (input.types.size() == 0)
    {
        std::cerr << "Error: Can't POSCAR file with no atoms.\n";
        return false;
    }

    std::ofstream outfile(filename);
    if (!outfile)
    {
        std::cerr << "Failed opening file : " << filename << " for read."
                  << std::endl;
        return false;
    }

    outfile.precision(std::numeric_limits<double>::digits10);
    outfile << "POSCAR generated by sp2\n"     // comment line
            << "1.0\n"; // lattice scaling

    // lattice vectors
    for (const auto &row : input.lattice)
    {
        for (double coord : row)
            outfile << coord << ' ';
        outfile << "\n";
    }

    auto temp = input;
    sort_structure_types(temp);

    // output types and figure out how many atoms of each type there are
    outfile << enum_to_str(temp.types[0]);
    std::vector<std::size_t> type_numbers = {1};
    for (std::size_t i = 1; i < temp.types.size(); ++i)
    {
        if (temp.types[i - 1] != temp.types[i])
        {
            outfile << ' ' << enum_to_str(temp.types[i]);
            type_numbers.push_back(1);
        }
        else
            type_numbers.back()++;
    }
    outfile << '\n';

    // output type numbers
    for (auto n : type_numbers)
        outfile << n << ' ';
    outfile << '\n';

    // x = L^T x'
    // (L^T)^-1 x = x'
    auto inv_lattice = mat3x3_t(input.lattice).inverse().transposed();
    outfile << "Direct\n";

    for (const auto &pos : temp.positions)
    {
        for (auto x : inv_lattice * pos)
            outfile << x << ' ';
        outfile << '\n';
    }

    return false;
}
