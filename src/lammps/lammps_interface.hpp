#ifndef SP2_LAMMPS_INTERFACE_HPP
#define SP2_LAMMPS_INTERFACE_HPP

#ifndef LAMMPS_INTERFACE_HPP_INCLUDED
#define LAMMPS_INTERFACE_HPP_INCLUDED

/// \file lammps_interface.hpp
/// \brief Interface header for the LAMMPS \cite plimpton1995fast library

#include "common/util/mpi.hpp"
#include "lammps/settings_t.hpp"
#include "common/enums.hpp"
#include "atac/atac_subsys_t.hpp"

#include <mpi.h>

#include <vector>
#include <string>

// forward declare of LAMMPS, to avoid having to include files here
/// LAMMPS \cite plimpton1995fast namespace
namespace LAMMPS_NS {
/// primary class for LAMMPS \cite plimpton1995fast
class LAMMPS;
}

namespace sp2 {
/// LAMMPS \cite plimpton1995fast interface namespace
namespace lammps {

/// Primary interface class to LAMMPS \cite plimpton1995fast
class system_control_t : public atac::atac_subsys_t
{
private:
    int na, ///< number of atoms
        nb; ///< number of bonds
    double total_potential;       ///< total system potential as calculated by LAMMPS
    std::vector<atom_type> type;  ///< atom type vector, defined in util_enums.hpp
    std::vector<double> position, ///< atom position vector
        force;                    ///< atom force vector

    /// lattice vectors
    double lattice[3][3];

    /// the only interface into LAMMPS, a pointer to a constructed LAMMPS class
    LAMMPS_NS::LAMMPS *lmp;

    /// MPI communicator information
    util::mpi_group_t mpi_info;

public:
    /// constructor, pulls communicator from MPI_COMM_WORLD, constructs lmp object
    /// \param mpi_group util::mpi_group_t input MPI group information and comms
    system_control_t(util::mpi_group_t mpi_group);
    /// destructor, automatically destroys the lammps object pointed to by lmp
    ~system_control_t();

    /// initialize the system
    void init(const structure_t &info,
        const lammps_settings_t &lmp_set);

    /// update forces and the total potential, executes a single LAMMPS run
    void update();
    /// pass a single string command to LAMMPS
    /// \param cmd string input command
    void exec_command(std::string cmd);

    template<class ...Args>
    void exec_command(std::string current, Args... other)
    {
        exec_command(current);
        exec_command(other...);
    }

    /// execute a command requiring a group, note that the formation of the group is temporary
    /// \param cmd string command to execute, the location of the group name in the command should be marked with %s
    /// \param ids const vector<int>& ids of the atoms that the command will be applied to
    void exec_group_command(std::string cmd, const std::vector<int> &ids);

    /// get the system potential
    double get_value() const;
    /// get the size of the gradient/position vectors
    size_t size() const;

    /// get atom positions
    std::vector<double> get_position() const;
    /// get atom forces
    std::vector<double> get_gradient() const;
    /// set atom positions
    /// \param input const vector<double>& input atom positions
    void set_position(const std::vector<double> &input);

    /// write atom positions/types to an xyz file
    /// \param filename string filename
    void write_output(std::string filename);
    /// append atom positions/types to an xyz file
    /// \param filename string filename
    void append_output(std::string filename);

    /// remove an atom
    void remove_atom(int id);

    /// add an atom
    void add_atom(atom_type type_in, const double *pos);

    /// conversion to bool for determining if everything is OK
    operator bool() {return lmp != nullptr;}

    /// conversion to structure interchange class
    structure_t get_structure() const;
    void set_structure(const structure_t  &input);

private:
    /// disallow copy constructor
    system_control_t(const system_control_t &);

    /// update number of atoms, positions, and types
    void update_atoms();

    /// update periodic boundary conditions
    void update_boundaries(bool initial, bool fixed = false);
};

} // namespace lammps
} // namespace sp2

#endif // SP2_LAMMPS_INTERFACE_HPP
