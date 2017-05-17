//
// Created by cc on 9/20/16.
//

#include "common/io/file_types/poscar.hpp"
#include "common/io/util.hpp"
#include "common/neighbor/utility_functions.hpp"

#include <string>
#include <fstream>
#include <common/math/vec3_t.hpp>

using namespace std;
using namespace sp2;

bool io::read_poscar(std::string filename, structure_t &output)
{
    try {
        ifstream infile(filename);
        if (!infile)
            throw runtime_error("could not open file");

        string line;

        // first line is a comment
        if (!getline(infile, line))
            throw runtime_error("could not read first line");

        // second line is lattice scaling factor
        if (!getline(infile, line))
            throw runtime_error("could not read second line");

        double lattice_scale;
        try {
            lattice_scale = stod(line);
        } catch (std::exception) {
            throw runtime_error("could not parse lattice scaling factor");
        }

        for (int i = 0; i < 3; ++i)
        {
            if (!getline(infile, line))
                throw runtime_error("could not read lattice vectors");

            auto split = io::split_line(line);
            if (split.size() != 3)
                throw runtime_error("could not parse lattice vectors");

            try {
                for (int j = 0; j < 3; ++j)
                    output.lattice[i][j] = lattice_scale * stod(split[j]);
            } catch (std::exception) {
                throw runtime_error("could not parse lattice vectors");
            }
        }

        if (!getline(infile, line) || line.empty())
            throw runtime_error("could not read atom types/numbers");

        // see if the next line is the number of atoms, or atom types
        bool has_atom_types = all_of(line.begin(), line.end(), [](char c) {
            return isspace(c) || isalpha(c);
        });

        vector<atom_type> types;
        if (!has_atom_types)
            types = {atom_type::CARBON};
        else
        {
            for (string type_str : io::split_line(line))
                types.push_back(enum_from_str<atom_type>(type_str));

            if (!getline(infile, line) || line.empty())
                throw runtime_error("could not read atom numbers");;
        }

        int total_atoms = 0;
        vector<int> n_atoms = {0};
        try {
            io::trim(line, " \t");
            for (string num_str : io::split_line(line))
                n_atoms.push_back(n_atoms.back() + stoi(num_str));

            total_atoms = n_atoms.back();
        } catch (std::exception) {
            throw runtime_error("could not parse atom numbers");
        }

        if (n_atoms.size() == 3 && types.size() == 1)
            types = {atom_type::HYDROGEN, atom_type::CARBON};

        if (!getline(infile, line) || line.empty())
            throw runtime_error("could not read coordinate specifier "
                "(direct/cartesian) or selective dynamics specifier");

        bool selective_dyn;
        if (std::toupper(line[0]) == 'S')
        {
            selective_dyn = true;
            if (!getline(infile, line) || line.empty())
                throw runtime_error("could not read coordinate specifier "
                    "(direct/cartesian)");
        }
        else
            selective_dyn = false;

        bool direct = std::toupper(line[0]) != 'C' &&
                      std::toupper(line[0]) != 'K';

        int type_n = 0;
        for (int n = 0; n < total_atoms && getline(infile, line); ++n)
        {
            auto positions = io::split_line(line);

            vec3_t pos = {};
            for (int i = 0; i < 3; ++i)
            {
                double val = stod(positions[i]);
                if (!direct)
                    pos[i] = val;
                else
                {
                    for (int j = 0; j < 3; ++j)
                        pos[j] += val * output.lattice[i][j];
                }
            }

            output.positions.push_back(pos);
            if (n >= n_atoms[type_n + 1])
                type_n++;

            output.types.push_back(types[type_n]);
        }
    } catch (const std::runtime_error ex) {
        cout << "Error: " << ex.what()
             << " when reading POSCAR file " << filename << endl;

        return false;
    } catch (const std::invalid_argument) {
        cout << "Error, failed parsing value when reading POSCAR file "
             << filename << endl;

        return false;
    }

    return true;
}

bool io::write_poscar(std::string filename, const structure_t &input)
{
    ofstream outfile(filename);
    if (!outfile)
    {
        cout << "Failed opening file : " << filename << " for read." << endl;
        return false;
    }

    outfile.precision(15);
    outfile << "\n"     // comment line
            << "1.0\n"; // lattice scaling

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
            outfile << input.lattice[i][j] << ' ';
        outfile << "\n";
    }

    auto temp = input;
    sort_structure_types(temp);

    int n_hydrogen = 0;
    for (auto t : temp.types)
        if (t != atom_type::CARBON)
            ++n_hydrogen;

    int n_carbon = temp.types.size() - n_hydrogen;
    auto positions = temp.positions;

    if (!n_hydrogen)
        outfile << "C\n" << n_carbon << '\n';
    else
        outfile << "H C\n" << n_hydrogen << ' ' << n_carbon << '\n';

    double inv_lattice[3][3];

    // x = L^T x'
    // (L^T)^-1 x = x'
    fbc::invert_3x3(input.lattice, inv_lattice);
    outfile << "Direct\n";

    for (std::size_t i = 0; i < positions.size(); ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            double val = 0;
            for (int k = 0; k < 3; ++k)
                val += positions[i][k] * inv_lattice[k][j];

            outfile << val << ' ';
        }
        outfile << '\n';
    }

    return false;
}
