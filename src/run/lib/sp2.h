#ifndef SP2_HEADER_H_INCLUDED
#define SP2_HEADER_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif

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
);

// same as calc_potential, except relaxes the entire structure + outputs
// the relaxed positions along with the other data
void relax_structure(
    int n_atoms,
    double *positions, // input/output in Angstroms
    double *gradient,  // output in eV / Angstrom
    double *potential, // output in eV
    double lattice[3][3]
);

#ifdef __cplusplus
}
#endif

#endif // SP2_HEADER_H_INCLUDED
