#include "common/io/file_types/cif.hpp"
#include "common/io/util.hpp"

#include <iostream>
#include <fstream>

using namespace std;
using namespace sp2;

bool io::write_cif(std::string filename, const structure_t &structure)
{
    const string info = "io::write_cif(): filename: \"" + filename + "\"";

    ofstream outfile(filename);
    if (!outfile.is_open())
    {
        cout << "Error opening file. " << info << endl;
        return false;
    }

    outfile.precision(8);
    outfile << "data_global\n"
            << "_chemical_name 'name'\n"
            << "_cell_length_a " << structure.lattice[0][0] << '\n'
            << "_cell_length_b " << structure.lattice[1][1] << '\n'
            << "_cell_length_c " << structure.lattice[2][2] << '\n'
            << "_cell_angle_alpha 90\n"
            << "_cell_angle_beta  90\n"
            << "_cell_angle_gamma 90\n"
            << "_symmetry_space_group_name_H-M '"
            << structure.space_group << "'\n"
            << "loop_\n"
            << "_atom_site_label\n"
            << "_atom_site_fract_x\n"
            << "_atom_site_fract_y\n"
            << "_atom_site_fract_z\n";

    auto n_atoms = structure.types.size();
    if (structure.n_symm > 0)
        n_atoms /= structure.n_symm;

    for (size_t i = 0; i < n_atoms; ++i)
    {
        outfile << enum_to_str(structure.types[i]);
        for (size_t j = 0; j < 3; ++j)
            outfile << '\t' << structure.positions[i * 3 + j]
                               / structure.lattice[j][j];
        outfile << '\n';
    }

    return true;
}

