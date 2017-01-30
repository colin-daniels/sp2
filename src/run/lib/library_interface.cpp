#include "sp2.h"

#include "common/io/util.hpp"
#include "common/io/structure.hpp"
#include "airebo/system_control_t.hpp"
#include "common/util/blas.hpp"
#include "common/minimize/minimize.hpp"
#include "common/vec3_t.hpp"
#include "common/json/json.hpp"

using namespace sp2;

structure_t get_structure(int n_atoms, double *positions, double lattice[3][3])
{
    structure_t structure;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            structure.lattice[i][j] = lattice[i][j];

    structure.types = std::vector<atom_type>(n_atoms, atom_type::CARBON);
    structure.positions = std::vector<double>(
        positions, positions + n_atoms * 3);

    return structure;
}

// calculate 2nd gen REBO potential
void calc_potential(
    int n_atoms,

    // For positions[] and gradient[]:
    // length: (n_atoms * 3)
    // layout: {x1, y1, z1, x2, y2, z2, ... }
    //
    // Calling code is responsible for allocation/lifetime and correct size.
    double *positions, // input in Angstroms
    double *gradient,  // output in eV / Angstrom
    double *potential, // output in eV

    // rows are lattice vectors
    double lattice[3][3]
) {
    airebo::system_control_t sys;
    sys.init(get_structure(n_atoms, positions, lattice));

    // copy gradient to output
    std::copy_n(
        sys.get_gradient().begin(),
        n_atoms * 3,
        gradient
    );

    // output potential
    *potential = sys.get_value();
}

void relax_structure(
    int n_atoms,
    double *positions, // input/output in Angstroms
    double *gradient,  // output in eV / Angstrom
    double *potential, // output in eV
    double lattice[3][3]
) {
    airebo::system_control_t sys;
    sys.init(get_structure(n_atoms, positions, lattice));

    minimize::acgsd_settings_t settings;
    settings.output_level = 0;
    settings.gradient_tolerance = 1e-4;
    minimize::acgsd(sys.get_diff_fn(), sys.get_position(), settings);

    // output gradient, position, and potential
    std::copy_n(sys.get_gradient().begin(), n_atoms * 3, gradient);
    std::copy_n(sys.get_position().begin(), n_atoms * 3, positions);
    *potential = sys.get_value();
}