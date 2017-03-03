//
// Created by cc on 9/9/16.
//

#ifndef SP2_AIREBO_SYSTEM_CONTROL_T_HPP
#define SP2_AIREBO_SYSTEM_CONTROL_T_HPP


/// \file system_control_t.hpp
/// \brief Main header for the REBO \cite brenner2002second \cite stuart2000reactive system class

#include "common/structure_t.hpp"
#include "common/neighbor/bond_control_t.hpp"
#include "common/math/blas.hpp"
#include "common/function_types.hpp"

#include <vector>
#include <string>

namespace sp2 {

/// REBO \cite brenner2002second \cite stuart2000reactive namespace
namespace airebo {

/// Main system class for the REBO \cite brenner2002second \cite stuart2000reactive potential
class system_control_t
{
private:
    int na,                         ///< number of atoms
        nb;                         ///< number of bonds
    double total_potential;         ///< total system potential energy
    double ref_lattice[3][3];       ///< original input lattice

    // atom info
    std::vector<atom_type> types;   ///< atom type vector (na x 1)
    std::vector<double> position,   ///< atom position vector (na x 3)
        force;                      ///< atom force vector (na x 3)

    std::vector<double> potential;  ///< per-atom potential

    std::vector<bool> static_atom;  ///< static atom flag vector, designates which atoms do not move/have 0 force (na x 1)
    std::vector<double> ref_pos;    ///< reference position vector used with the static atom vector (na x 3)

    /// bond controller class, used to update connectivity and calculate basic bond information
    fbc::bond_control_t bond_control;

    std::vector<bond_type> b_types; ///< bond type vector (e.g. carbon-hydrogen, carbon-carbon, etc)
    std::vector<double> bond_force, ///< bond force vector (nb x 3)
        bond_dir_force,  ///< vector for forces for the individual bond directions (nb x 1)
        force_on_cutoff; ///< forces that multiply the cutoff function of a given bond (nb x 1)

    std::vector<double> T,          ///< temporary storage for values of T (nb x 4), (T, dT/dN_ti, dT/dN_tj, dT/dN_conj)
        F,                          ///< temporary storage for values of F (nb x 4), (F, dF/dN_ti, dF/dN_tj, dF/dN_conj)
        P;                          ///< temporary storage for values of P (nb x 3), (P, dP/dN_h, dP/dN_c)
    std::vector<double> Na_vars,    ///< atom information vector for number of bonds (na x 3), (N_t, N_h, N_c, N_conj_half)
        Nb_vars;                    ///< bond information for connected atom's number of bonds (nb x 7), (N_ti, N_tj, N_h, N_c, N_conj, dNca/dr_ij, dNcb/dr_ij)

    std::vector<double> cutoff,     ///< cutoff function values (nb x 2), (cutoff, d/dr cutoff)
        attr_p;                     ///< attractive potential values (nb x 2), (attr_p, d/dr attr_p)

    std::vector<double> self_dots,  ///< inner-atom dot product storage vector
        other_dots,                 ///< inter-atom dot product storage vector
        cos_vals;                   ///< cosine values

    std::vector<int> sd_offsets,    ///< offsets for inner-atom dot products
        od_offsets;                 ///< offsets for inter-atom dot products

    bool vdw_enabled = false;
    std::vector<double> vdw_deltas, ///< temp deltas
        vdw_force,                  ///< temporary storage for van der Waals force values, to avoid re-allocating on every update
        vdw_working;                ///< working storage for van der Waals calculations

public:
    system_control_t() = default;
    system_control_t(const system_control_t &) = delete;
    system_control_t& operator=(const system_control_t&) = delete;

    /// \brief initialization function
    /// \param lattice double[3][3] input lattice vector matrix (rows) to designate periodicity, rows with all elements zero are assumed non-periodic
    /// \param position_in std::vector<double> input positions
    /// \param types_in std::vector<atom_type> input types
    void init(const double lattice[3][3],
        const std::vector<double> &position_in,
        const std::vector<atom_type> &types_in);

    /// \brief initialization function
    void init(const structure_t &structure);

    /// set/update lattice vectors
    void set_lattice(const double lattice[3][3]);

    /// wrap all atom positions into the current unit cell
    void wrap_positions();

    /// update function, calculates total potential and forces for the current atom positions
    void update();

    /// get the number of atoms in the system
    int get_na() const {return na;}
    /// get the number of bonds in the system
    int get_nb() const {return nb;}
    /// get the number of active cells in the system's bond controller class
    int get_nc() const {return bond_control.get_nc();}

    /// get the current lattice
    void get_lattice(double output[3][3]) {
        std::copy_n(ref_lattice[0], 9, output[0]);}
    /// get a reference to the system's bond controller class
    fbc::bond_control_t &get_bond_control() {return bond_control;}

    /// get the position vector
    std::vector<double> get_position() const {return position;}
    /// get the force vector
    std::vector<double> get_gradient() const {
        std::vector<double> grad(force); vscal(-1.0, grad); return grad;}
    /// get the atom type vector
    std::vector<atom_type> get_types() const {return types;}
    /// get the bond force vector
    std::vector<double> get_bond_forces() const {return bond_force;}
    /// get the static atom vector
    std::vector<bool> get_static_atoms() const {return static_atom;}
    /// get per-atom potentials
    std::vector<double> get_potentials() const {return potential;}

    /// set the position vector, does not update the reference position vector used with static atoms
    void set_position(const std::vector<double> &input_pos) {
        position = input_pos;}
    /// set the type vector
    void set_types(const std::vector<atom_type> &types_in) {
        types = types_in;}
    /// designate which atoms are static in the system
    void set_static_atoms(const std::vector<bool> &static_atoms_in) {
        static_atom = static_atoms_in;}
    /// update the reference position used for static atoms, effectively moves any static atoms
    void set_position_static(const std::vector<double> &input_pos) {
        position = input_pos; ref_pos = input_pos;}

    /// get the total size of the system, equal to the length of the force or position vectors
    size_t size() const {return position.size();}
    /// get the total potential
    double get_value() const {return total_potential;}
    /// get the pointer to the beginning of the force vector
    double *force_ptr() {return force.data();}
    /// get the pointer to the beginning of the position vector
    double *position_ptr() {return position.data();}

    /// write atom positions and types to a file
    void write_output(std::string filename);
    /// append atom positions and types to a file
    void append_output(std::string filename);

    /// \brief add hydrogen to carbon atoms
    void add_hydrogen();

    structure_t get_structure() const;
    void set_structure(const structure_t &input);

    /// construct a differentiable function object tied to this object
    diff_fn_t get_diff_fn();

private:
    /// update connectivity and bond lengths/components
    void bond_update();
    /// calculate cutoff functions, attractive/repulsive potential terms, and fill Na_vars
    void update_stage1();
    /// calculate entries for Nb_vars
    void update_stage2();
    /// calculate T, P, and F vectors
    void update_stage3();
    /// calculate the inner and inter-atom dot products
    void update_stage4();
    /// calculate G(theta) terms and b^{sigma - pi}_ij
    void update_stage5();
    /// calculate torsion terms in b^pi_ij
    void update_stage6();
    /// calculate derivatives for N_conj related quantities
    void update_stage7();
    /// distribute calculated forces from bonds to atoms
    void update_distribute();

    void update_vdw();

    /// TODO: calculate forces on the lattice vectors
    void update_lattice_force();
};

} // namespace rebo
} // namespace sp2

#endif // SP2_AIREBO_SYSTEM_CONTROL_T_HPP
