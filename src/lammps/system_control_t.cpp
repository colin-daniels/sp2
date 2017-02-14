#include "lammps_interface.hpp"
#include "common/io/util.hpp"
#include "common/io/structure.hpp"
#include "common/util/templates.hpp"
#include "common/structure_t.hpp"

#include <sstream>

// lammps include files
#include <lammps/lammps.h>
#include <lammps/input.h>
#include <lammps/atom.h>
#include <lammps/library.h>

#include <boost/mpi/communicator.hpp>
#include <common/util/blas.hpp>
#include <common/io/util.hpp>
#include <common/vec3_t.hpp>
#include "common/neighbor/utility_functions.hpp"

using namespace LAMMPS_NS;

using namespace std;
using namespace sp2;

void transform_input(std::vector<double> &pos, double lattice[3][3],
    double inv_rotation[3][3])
{
    vec3_t lat[3] = {
        vec3_t(lattice[0]),
        vec3_t(lattice[1]),
        vec3_t(lattice[2])
    };

    if (dot(cross(lat[0], lat[1]), lat[2]) < 0)
    {
        // not a right-handed basis
        // swap rows to fix
        std::swap(lat[0], lat[1]);
        for (int i = 0; i < 3; ++i)
            std::swap(lattice[0][i], lattice[1][i]);
    }

    double a = lat[0].mag(),
        b = lat[1].mag(),
        c = lat[2].mag();

    for (auto &v : lat)
        v.normalize();

    double cos_alpha = dot(lat[1], lat[2]), // b/c
        cos_beta = dot(lat[0], lat[2]),     // a/c
        cos_gamma = dot(lat[0], lat[1]);    // a/b

    // http://lammps.sandia.gov/doc/Section_howto.html#howto-12
    double lx = a,
        xy = b * cos_gamma,
        xz = c * cos_beta,
        ly = sqrt(b * b - xy * xy),
        yz = (b * c * cos_alpha - xy * xz) / ly,
        lz = sqrt(c * c - xz * xz - yz * yz);

    double new_lattice[3][3] = {
        {lx,  0,  0},
        {xy, ly,  0},
        {xz, yz, lz}
    };

    // get rotation matrix R via solving
    //     L R = N
    double inv_lattice[3][3] = {};
    fbc::invert_3x3(lattice, inv_lattice);

    // rotation matrix should be its transpose, because we are going to use it
    // to rotate all of the vectors next
    double rotation[3][3] = {};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                rotation[j][i] += inv_lattice[i][k] * new_lattice[k][j];

    // copy lattice to output
    std::copy_n(new_lattice[0], 9, lattice[0]);

    // copy transpose rotation to inv_rotation
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            inv_rotation[i][j] = rotation[j][i];

    // transform positions
    for (std::size_t i = 0; i < pos.size() / 3; ++i)
    {
        sp2::vec3_t new_pos = vec3_t(&pos[i * 3]).mul_3x3(rotation);
        for (int j = 0; j < 3; ++j)
            pos[i * 3 + j] = new_pos[j];
    }
}

lammps::system_control_t::system_control_t(boost::mpi::communicator comm) :
    na(0), nb(0), total_potential(0), lmp(nullptr),
    comm(comm) {}

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
    sp2::vscal(-1.0, grad);
    return grad;
}

void lammps::system_control_t::set_position(const vector<double> &input) {
    position = input;}

void lammps::system_control_t::set_structure(const structure_t  &input)
{
    type = input.types;
    position = input.positions;

    bool eq = true;
    for (int i = 0; i < 9; ++i)
        eq &= input.lattice[i%3][i/3] == lattice_orig[i%3][i/3];

    if (!eq)
    {
        std::copy_n(input.lattice[0], 9, lattice_orig[0]);
        transform_input(position, lattice, inv_transform);
        update_boundaries(false);
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

void lammps::system_control_t::init(const structure_t &info,
    const lammps_settings_t &lmp_set)
{
    copy_n(info.lattice[0], 9, lattice[0]);

    // number of atoms
    na = info.types.size();

    type = info.types;
    position = info.positions;
    force.resize(na * 3, 0);

    // transform to triclinic
    std::copy_n(lattice[0], 9, lattice_orig[0]);
    transform_input(position, lattice, inv_transform);

    delete lmp;

    if (lmp_set.log_output)
        lmp = new LAMMPS(0, NULL, static_cast<MPI_Comm>(comm));
    else
    {
        // LAMMPS constructor takes char* not const char*
        char dir[] = "/",
            opt1[] = "-screen",
            opt2[] = "none";
        char *args[] = {dir, opt1, opt2};
        lmp = new LAMMPS(3, args, static_cast<MPI_Comm>(comm));
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

    comm.barrier();

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
void lammps::system_control_t::write_output(std::string filename)
{
    if (comm.rank() == 0)
        io::write_structure(filename, get_structure());
}

void lammps::system_control_t::append_output(std::string filename)
{
    if (comm.rank() == 0)
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
        // non-periodic?
        // TODO: fix for triclinic
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

#warning TEMP

    // tell lammps what directions are periodic
    exec_command(periodicity_cmd);

    std::ostringstream ss;
    for (int i = 0; i < 3; ++i)
        ss << "0 " << lattice[i][i] << ' ';

    // if triclinic
    ss << lattice[1][0] << ' ' << lattice[2][0] << ' ' << lattice[2][1];

    if (initial)
    {
        exec_command(
            "box tilt large",
            // define the extents of the simulation box
            "region sim prism " + ss.str(),
            // create the simulation box
            "create_box 2 sim" // note the 2 is because there are 2 atom types
        );
    }
    else
    {
//        // modify simulation box extents
//        exec_command(
//            "change_box all x final "
//            + to_string(min_max[0][0]) + " "
//            + to_string(min_max[0][1]),
//            "change_box all y final "
//            + to_string(min_max[1][0]) + " "
//            + to_string(min_max[1][1]),
//            "change_box all z final "
//            + to_string(min_max[2][0]) + " "
//            + to_string(min_max[2][1])
//        );
    }
}

#ifdef SP2_DEBUG

#include <gtest/gtest.h>
#include "common/io/structure.hpp"

TEST(lammps, all) {
    sp2::structure_t input;
    io::read_structure("graphene.xyz", input);

    lammps::system_control_t sys;
    sys.init(input);

    sys.update();

    EXPECT_NE(sys.get_value(), 0.0);
}

#endif // SP2_DEBUG
