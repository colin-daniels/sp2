#include "lammps_interface.hpp"
#include "common/io/util.hpp
#include "common/io/structure.hpp"
#include "common/util/templates.hpp"
#include "common/structure_t.hpp"

#include <fstream>

// lammps include files
#include <lammps/lammps.h>
#include <lammps/input.h>
#include <lammps/atom.h>
#include <lammps/library.h>

#include <cassert>

using namespace LAMMPS_NS;

using namespace std;
using namespace sp2;

lammps::system_control_t::system_control_t(util::mpi_group_t mpi_group) :
    na(0), nb(0), total_potential(0), lmp(nullptr),
    mpi_info(mpi_group) {}

lammps::system_control_t::~system_control_t() {
    delete lmp;}

double lammps::system_control_t::get_value() const {
    return total_potential;}

size_t lammps::system_control_t::size() const {
    return position.size();}

vector<double> lammps::system_control_t::get_position() const {
    return position;}

vector<double> lammps::system_control_t::get_gradient() const
{
    vector<double> grad(force);
    vscal(-1, grad);
    return grad;
}

void lammps::system_control_t::set_position(const vector<double> &input) {
    position = input;}

void lammps::system_control_t::set_structure(const structure_t  &input)
{
    type = input.types;
    position = input.positions;
    for (int i = 0; i < 9; ++i)
    {
        if (input.lattice[i/3][i%3] == lattice[i/3][i%3])
            continue;

        copy_n(input.lattice[0], 9, lattice[0]);
        update_boundaries(false);
        break;
    }
}

void lammps::system_control_t::init(const structure_t &info,
    const lammps_settings_t &lmp_set)
{
    copy_n(info.lattice[0], 9, lattice[0]);

    // number of atoms
    na = info.types.size();

    type = info.types;
    position = info.positions;
    force.resize(na * 3, 0);

    if (!io::file_exists("CH.airebo"))
    {
        cerr << "CH.airebo is missing, cannot use LAMMPS, aborting..." << endl;
        abort();
    }

    delete lmp;

    if (lmp_set.log_output)
        lmp = new LAMMPS(0, NULL, mpi_info.group_comm);
    else
    {
        // LAMMPS constructor takes char* not const char*
        char dir[] = "/",
            opt1[] = "-screen",
            opt2[] = "none";
        char *args[] = {dir, opt1, opt2};
        lmp = new LAMMPS(3, args, mpi_info.group_comm);
    }

    if (!lmp)
        return;

    exec_command(
        "units metal",                  // Angstroms, picoseconds, eV
        "processors * * *",             // automatic processor mapping
        "atom_style atomic",            // attributes to store per-atom
        "thermo_modify lost error",     // don't let atoms disappear
        // without telling us
        "atom_modify map array",        // store all positions in an array
        "atom_modify sort 0 0.0"        // don't reorder atoms during simulation
    );

    // set up simulation box and periodicity
    update_boundaries(true /* initial */);

    // set atom masses and add them to the system
    exec_command("mass 1 1.0", "mass 2 12.01");
    for (int i = 0; i < na; ++i)
    {
        exec_command("create_atoms "
                     + to_string((int)type[i] + 1) + " single" // lammps types start at 1
                     + " " + to_string(position[i * 3])
                     + " " + to_string(position[i * 3 + 1])
                     + " " + to_string(position[i * 3 + 2])
                     + " remap yes" // map periodic coordinates back into the box
        );
    }

    MPI_Barrier(mpi_info.group_comm);

    exec_command(
        "pair_style airebo "            // airebo pair style potential
        + to_string(lmp_set.sigma_scale) + " "          // LJ range (x3.4 A)
        + string(lmp_set.compute_lj ? "1 " : "0 ")      // LJ on/off
        + string(lmp_set.compute_torsion ? "1" : "0"),  // torsion on/off
        "pair_coeff * * CH.airebo H C", // read potential info
        "compute 1 all pe",             // set up compute ID 1 for energy
        "run 0"                         // do an iteration
    );

    assert(static_cast<int>(lmp->atom->natoms) == na);
}

void lammps::system_control_t::update()
{
    // update lammps with any changes to number of atoms, positions, or types
    update_atoms();

    // run lammps for one interation (run 1 is two iterations)
    exec_command("run 0");

    // get the force
    char cf[] = "f";
    lammps_gather_atoms(lmp, cf, 1, 3, force.data());

    // and the potential (defined in init() as the first compute)
    char comp_n[] = "1";
    total_potential = *static_cast<double*>(
        lammps_extract_compute(lmp, comp_n, 0, 0));
}

void lammps::system_control_t::update_atoms()
{
    // get the number of atoms from lammps
    std::size_t lmp_na = static_cast<std::size_t>(lmp->atom->natoms);

    // update positions
    char cx[] = "x";
    lammps_scatter_atoms(lmp, cx, 1, 3, position.data());

    // change any different atom types individually
    // since I'm not sure if scattering types has any side effects
    char ct[] = "type";
    vector<int> lmp_type(lmp_na, 0);
    lammps_gather_atoms(lmp, ct, 0, 1, lmp_type.data());

    for (std::size_t i = 0; i < min(lmp_type.size(), type.size()); ++i)
        if (lmp_type[i] != (int)type[i] + 1) // (type-ids start at 1)
            exec_command("set atom " + to_string(i)
                         + " type " + to_string((int)type[i] + 1));

    // if atoms were removed since the last update
    if (type.size() < lmp_na)
    {
        vector<int> ids_to_delete;
        for (std::size_t i = type.size(); i < lmp_na; ++i)
            ids_to_delete.push_back(i);

        exec_group_command("delete_atoms group %s", ids_to_delete);
    }
    else if (type.size() > lmp_na)
    {
        // if atoms were added since the last update, add them
        update_boundaries(false /* not initial */, true  /* fixed */);

        for (std::size_t i = lmp_na; i < type.size(); ++i)
            exec_command("create_atoms "
                         + to_string((int)type[i] + 1) + " single"
                         + " " + to_string(position[i * 3])
                         + " " + to_string(position[i * 3 + 1])
                         + " " + to_string(position[i * 3 + 2])
                         + " remap yes");

        // reset non-periodic bounds to shrink-wrap mode
        update_boundaries(false /* not initial */, false /* shrink-wrap */);
    }
}

structure_t lammps::system_control_t::get_structure() const
{
    structure_t structure;

    structure.types = type;
    structure.positions = position;
    copy_n(lattice[0], 9, structure.lattice[0]);

    return structure;
}

void lammps::system_control_t::write_output(std::string filename)
{
    if (mpi_info.group_rank == 0)
        io::write_structure(filename, get_structure());
}

void lammps::system_control_t::append_output(std::string filename)
{
    if (mpi_info.group_rank == 0)
        io::write_structure(filename, get_structure(), true);
}

void lammps::system_control_t::add_atom(atom_type type_in, const double *pos)
{
    stringstream ss;
    ss << "create_atoms " << (int)type_in + 1 << " single";
    for (int i = 0; i < 3; ++i)
    {
        // write the addition command
        ss << ' ' << pos[i];

        // resize vectors/record input
        position.push_back(pos[i]);
        force.push_back(0);
    }

    // increase number of atoms
    na += 1;
    type.push_back(type_in);

    // need to update periodic bounds (and set to be fixed first)
    update_boundaries(false /* not initial */, true  /* fixed bounds */);

    // add the atom, making sure to wrap periodic coordinates (remap option)
    exec_command(ss.str() + " remap yes");

    // reset to shrink wrap style
    update_boundaries(false /* not initial */, false /* shrink-wrap bounds */);
}


void lammps::system_control_t::remove_atom(int id)
{
    // swap atom with the last one
    for (int i = 0; i < 3; ++i)
    {
        swap(position[id * 3 + i], position[(na - 1) * 3 + i]);
        swap(force[id * 3 + i], force[(na - 1) * 3 + i]);
    }

    // execute the command to remove the last atom
    exec_group_command("delete_atoms group %s", {na - 1});

    // update number of atoms
    na -= 1;
    position.resize(na * 3);
    force.resize(na * 3);
}

void lammps::system_control_t::exec_command(std::string cmd)
{
    lmp->input->one(cmd.c_str());
}

void lammps::system_control_t::exec_group_command(std::string cmd,
    const std::vector<int> &ids)
{
    const string group_name = "tempg";
    if (ids.empty())
        return;

    // generate the group creation command
    string add_group = "group " + group_name + " id";
    for (size_t i = 0; i < ids.size(); ++i)
        add_group += " " + to_string(ids[i]);

    // make the temporary group
    exec_command(add_group);

    // insert the group name into the command
    size_t pos = cmd.find("%s");
    if (pos == string::npos)
        return;

    cmd.replace(pos, 2, group_name);

    // execute the command
    exec_command(cmd);

    // delete the temporary group
    exec_command("group " + group_name + " delete");
}


void lammps::system_control_t::update_boundaries(bool initial, bool fixed)
{
    // determine min/max positions of atoms in the system
    // and set up periodic boundary conditions
    using nld = numeric_limits<double>;
    double min_max[3][2] = {
        {nld::max(), nld::lowest()},
        {nld::max(), nld::lowest()},
        {nld::max(), nld::lowest()}
    };

    // need a bit extra for box padding because LAMMPS
    // doesn't like it if it's too small
    const double padding = 2;
    for (int i = 0; i < na; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            min_max[j][0] = std::min(min_max[j][0],
                position[i * 3 + j] - padding);
            min_max[j][1] = std::max(min_max[j][1],
                position[i * 3 + j] + padding);
        }
    }

    string periodicity_cmd;
    if (initial)
        periodicity_cmd = "boundary";
    else
        periodicity_cmd = "change_box all boundary";

    for (int i = 0; i < 3; ++i)
    {
        // non-periodic
        if (lattice[i][i] == 0)
            periodicity_cmd += (fixed ? " f" : " s");
        else
        {
            // perioidic
            periodicity_cmd += " p";

            // reset min/max to the lattice vector length
            min_max[i][0] = 0;
            min_max[i][1] = lattice[i][i];
        }
    }

    // tell lammps what directions are periodic
    exec_command(periodicity_cmd);

    if (initial)
    {
        exec_command(
            // define the extents of the simulation box
            "region sim block  "
            + to_string(min_max[0][0]) + " "
            + to_string(min_max[0][1]) + " "
            + to_string(min_max[1][0]) + " "
            + to_string(min_max[1][1]) + " "
            + to_string(min_max[2][0]) + " "
            + to_string(min_max[2][1]),
            // create the simulation box
            "create_box 2 sim" // note the 2 is because there are 2 atom types
        );
    }
    else
    {
        // modify simulation box extents
        exec_command(
            "change_box all x final "
            + to_string(min_max[0][0]) + " "
            + to_string(min_max[0][1]),
            "change_box all y final "
            + to_string(min_max[1][0]) + " "
            + to_string(min_max[1][1]),
            "change_box all z final "
            + to_string(min_max[2][0]) + " "
            + to_string(min_max[2][1])
        );
    }
}
