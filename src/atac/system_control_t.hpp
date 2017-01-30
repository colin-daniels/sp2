#ifndef SP2_ATAC_SYSTEM_CONTROL_T_HPP
#define SP2_ATAC_SYSTEM_CONTROL_T_HPP

#include "common/util/mpi.hpp"
#include "common/minimize/minimizeable_t.hpp"

#include "common/neighbor/bond_control_t.hpp"
#include "atac/settings_t.hpp"
#include "common/enums.hpp"

#include "atac/atac_subsys_t.hpp"

#include <memory>
#include <functional>
#include <vector>
#include <set>
#include <mpi.h>

namespace sp2 {
/// ATAC \cite bullard2013dynamical \cite daniels2014emergent library namespace
namespace atac {
/// Boltzmann constant in eV/K
static const double boltzmann = 8.6173e-5;

/// mutation structure type containing mutation name,
/// probability, and application function (should probably make
//  this a larger class)
struct mut_instance_t
{
    /// mutation abbreviation/name, must be unique
    std::string name;
    /// relative probability of application (not acceptance) [0 to 1]
    double probability;
    /// function which applies the mutation when called
    std::function<bool()> apply;

    /// whether to always accept the mutation regardless of the
    /// new system energy
    bool always_accept;

    /// comparison op for set storage
    bool operator<(const mut_instance_t &b) const {
        return name < b.name;}
};

/// primary class used with the ATAC algorithm
class system_control_t : public minimize::differentiable_t
{
private:
    int na, ///< number of atoms
        nb; ///< number of bonds

    int iteration,                  ///< current iteration number
        accepted_mutations;         ///< number of accepted mutations/iterations

    double best_delta;              ///< lowest mutation delta E since the last accepted mutation

    double total_potential,         ///< system potential
        penalty_potential;          ///< potential due to the penalty functions

    graph::ud_graph_t graph;        ///< graph representing enforced bonds
    std::vector<double> deltas,     ///< position data for penalty potentials
        lengths;                    ///< bond length data for penalty potentials

    double lattice[3][3];           ///< lattice vectors
    std::vector<atom_type> type;    ///< atom types
    std::vector<double> position,   ///< atom positions
        gradient;                   ///< atom gradient data

    util::mpi_group_t mpi_info;     ///< MPI communicator/group information
    atac_settings_t settings;      ///< general settings for ATAC

    /// set of possible mutations, populated in system_control_t::init()
    std::set<mut_instance_t> mutations;

    /// main input system, generates forces and potential
    std::unique_ptr<atac_subsys_t> sub_sys;

    /// bond controller, calculates adjacency information and
    /// bond vector components for penalty potential calculation
    fbc::bond_control_t bond_control;

public:
    system_control_t(util::mpi_group_t mpi_group);

    /// initialize the system and take ownership of sys_in
    void init(std::unique_ptr<atac_subsys_t> &&sys_in,
        atac_settings_t input_settings);

    /// execute one iteration of ATAC
    bool iterate();

    /// call the subsystem's update function and calculate penalty
    /// potential values and forces
    void update();

    /// get the total potential as calculated by the last call to update()
    double get_value() const {return total_potential;}
    /// get the size of the position/gradient vectors
    size_t size() const {return position.size();}

    /// get atom positions
    std::vector<double> get_position() const {return position;}
    /// get the gradient vector as calculated by the last call to update()
    std::vector<double> get_gradient() const {return gradient;}
    /// set the atom positions
    void set_position(const std::vector<double> &position_in) {
        position = position_in;}

    void set_structure(const structure_t &input);
    structure_t get_structure() const;

    void bcast(MPI_Comm comm, int root);

    void write_output(std::string filename) {
        sub_sys->write_output(filename);}
    void append_output(std::string filename) {
        sub_sys->append_output(filename);}

private:
    /// calculate energy and forces from the penalty potential
    void update_penalty();

    /// \brief repeatedly pick and attempt to apply a mutation to the system
    /// \returns type that was successfully applied (or an error type)
    mut_instance_t attempt_mutation(int max_attempts = 1000);

    /// randomly pick a type of mutation according to current settings
    /// note: called by system_control_t::attempt_mutation()
    mut_instance_t pick_mutation();

    /// mutation, Stone-Thrower-Wales bond rotation
    bool stone_wales();
    /// mutation, two atom vacancy (remove dimer)
    bool divacancy();
    /// mutation, add two atoms
    bool add_dimer();
    /// mutation, sp2->sp3->sp2 cycloaddition (two bonds)
    bool cycloaddition();

    bool stretch_lattice();

    /// remove a specific atom
    void remove_atom(int id);
    /// create an atom with a given type and position
    int create_atom(atom_type type_in, double *pos);

    /// basically every single mutation involves some kind of
    /// pair of sp2 bonded carbon atoms
    std::vector<graph::ud_edge_t> random_sp2_carbon_bonds();
};

} // namespace atac
} // namespace sp2

#endif // SP2_ATAC_SYSTEM_CONTROL_T_HPP
