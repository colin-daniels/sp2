#include "run/run_types.hpp"

#include "common/io/util.hpp"
#include "common/io/structure.hpp"

#include "common/minimize/minimize.hpp"
#include "common/minimize/pso/particle_t.hpp"

#include "symm/system_control_t.hpp"
#include "symm/space_group_t.hpp"
#include "symm/symm_settings_t.hpp"
#include "symm/pso_adapters.hpp"
#include "symm/util.hpp"

#include "common/neighbor/bond_control_t.hpp"

#include <iostream>
#include <boost/mpi.hpp>

using namespace std;
using namespace sp2;
namespace mpi = boost::mpi;

/// takes full structure, penalizes long bonds and neighbors != 3
double penalty_func(const structure_t &structure)
{
    /// penalty potential multiplier
    constexpr double multiplier = 0.0;

    static fbc::bond_control_t bond_control;
    bond_control.init(structure.lattice, 2.0, 0);

    bond_control.update(sp2::v3tod(structure.positions));

    auto graph = bond_control.get_graph();
    auto bond_lengths = bond_control.get_bond_lengths();

    double penalty_potential = 0;
    for (auto atom_id : graph.vertices())
    {
        double n_bonds = 0;
        for (auto bond : graph.edges(atom_id))
        {
            auto bond_len = bond_lengths[bond.id];
            if (bond_len < 1.55)
                n_bonds += 1;
            else if (bond_len < 1.65)
                n_bonds += (1.65 - bond_len) / 0.1;
        }

        auto diff = std::abs(n_bonds - 3.0);
        penalty_potential += diff * multiplier;
    }

    return penalty_potential;
}

void process_files()
{
    auto groups = symm::read_groups("space-groups.txt");
    cout.precision(10);

    symm::system_control_t sys;
    sys.set_group(groups["I a -3 d"]);

    minimize::acgsd_settings_t cg_set;
    cg_set.iteration_limit = 100;
    cg_set.value_tolerance = 1e-4;
    cg_set.output_level = 1;

    for (string filename : {
        "f_gyroid1a.xyz",
        "f_gyroid2a.xyz",
        "f_gyroid3a.xyz",
        "f_gyroid4a.xyz",
        "f_gyroid5a.xyz",
        "f_gyroid6a.xyz",
        "f_gyroid7a.xyz",
        "f_gyroid8a_0.xyz",
        "f_gyroid8a_1.xyz",
        "f_gyroid8a_2.xyz",
        "f_gyroid8a_3.xyz",
        "f_gyroid9a.xyz",
        "f_gyroid10a.xyz",
        "f_gyroid11a.xyz",
        "f_gyroid12a.xyz",
        "f_gyroid13a.xyz",
        "f_gyroid14a.xyz",
        "f_gyroid15a.xyz",
        "f_gyroid16a.xyz",
        "f_gyroid17a.xyz",
        "f_gyroid18a.xyz",
        "f_gyroid21a.xyz"
    }) {

        structure_t structure;
        if (!io::read_structure(filename, structure))
        {
            cerr << "Failed to read: " << filename << endl;
            continue;
        }
        auto uc = structure.lattice[0][0];
        sys.set_structure(structure);
        sys.update();

        auto ref_structure = structure;
        for (auto &v : ref_structure.positions)
            v /= uc;

        auto get_struct = [&](double a) {
            auto temp = ref_structure;
            for (auto &v : temp.positions)
                v *= a;

            for (int i = 0; i < 3; ++i)
                temp.lattice[i][i] = a;

            return temp;
        };

        auto get_value = [&](double a) {
            sys.set_structure(get_struct(a + uc));
            sys.update();

            if (vdot(sys.get_gradient(), sys.get_gradient()) < 1e-4)
                return sys.get_value();

            minimize::acgsd([&](const auto &pos) {
                auto temp = sys.get_structure();
                temp.positions = sp2::dtov3(pos);
                sys.set_structure(temp);

                sys.update();
                return make_pair(sys.get_value(), sys.get_gradient());
            }, sp2::v3tod(sys.get_structure().positions), cg_set);

            return sys.get_value();
        };

        auto get_slope = [&](double a) {
            constexpr double step = 1e-4;
            return (get_value(a + step) - get_value(a - step)) / (2 * step);
        };

        get_value(0);
        ref_structure.positions = sys.get_structure().positions;
        for (auto &v : ref_structure.positions)
            v /= uc;

        double alpha = get_slope(0);
        cout << filename << " v: " << sys.get_value()
             << "\tuc: " << uc << "\ts: " << alpha << endl;
        for (int i = 0; i < 5 && std::abs(alpha) > 1e-4; ++i)
        {
            alpha = minimize::linesearch(get_value, get_slope,
                max(min(alpha, 0.5), -0.5));

            if (alpha == 0)
                break;

            uc += alpha;
            alpha = get_slope(0);
            cout << filename << " v: " << sys.get_value() << "\tuc: "
                 << uc << "\ts: " <<  alpha << endl;
        }

        structure = sys.get_full_structure();
        for (auto &v : structure.positions)
            v -= uc * floor(v / uc);

        io::write_structure("c" + filename, structure);
        sys.update();
        cout << "grep:" << filename << " " << sys.get_value() << endl;
    }
}

int sp2::run_symm(const run_settings_t &settings, MPI_Comm comm_in)
{
    mpi::communicator comm(comm_in, mpi::comm_attach);

    ostream error_log(cerr.rdbuf());
    // only output if rank 0 (badbit prevents error_log from outputting)
    if (comm.rank() != 0)
        error_log.setstate(ios_base::badbit);

////////////////////////////////////////////////////////////////////////////////
// Load settings from input config file and do some sanity checks             //
////////////////////////////////////////////////////////////////////////////////

    symm::symm_settings_t symm_settings = settings.symm_settings;

    // number of atoms (irreducible representation)
    if (symm_settings.n_atoms <= 0)
    {
        error_log << "Error, n_atoms <= 0" << endl;
        return EXIT_FAILURE;
    }

    // unit cell size
    if (symm_settings.unit_cell_range[0] == symm_settings.unit_cell_range[1] &&
        symm_settings.unit_cell_range[0] == 0)
    {
        error_log << "Error, unit cell range equal to zero." << endl;
        return EXIT_FAILURE;
    }

////////////////////////////////////////////////////////////////////////////////
// Load space groups and find the specified one                               //
////////////////////////////////////////////////////////////////////////////////

    // load space groups
    if (!io::file_exists(symm_settings.space_group_file))
    {
        error_log << "Error, missing space group file \""
                  << symm_settings.space_group_file << "\"." << endl;
        return EXIT_FAILURE;
    }

    // check for the input name
    auto groups = symm::read_groups(symm_settings.space_group_file);
    if (!groups.count(symm_settings.space_group))
    {
        error_log << "Couldn't find space group with name \""
                  << symm_settings.space_group << "\" in space group file." << endl;
        return EXIT_FAILURE;
    }

////////////////////////////////////////////////////////////////////////////////
// Initialize potential calculation object (REBO)                             //
////////////////////////////////////////////////////////////////////////////////

    symm::system_control_t sys;
    cout << "space group: " << symm_settings.space_group << endl;
    sys.set_group(groups[symm_settings.space_group]);

    // function that will convert pso output to atom positions
    // (note: we could just do a search where the pso output is actual
    //  cartesian atom coordinates, but this is hopefully better)
    // auto pso_adapter = bind(symm::bonded_adapter,
    //     std::placeholders::_1,
    //     settings.unit_cell_range,
    //     settings.bond_range,
    //     settings.n_connected);

    auto pso_adapter = bind(symm::basic_adapter,
        std::placeholders::_1,
        symm_settings.unit_cell_range);

    auto inv_adapter = bind(symm::inverse_basic_adapter,
        std::placeholders::_1,
        symm_settings.unit_cell_range);

    // min/max bounds set to just be 1.0, we will generate the actual
    // atom positions from the pso output
    symm_settings.pso_set.min_bound = vector<double>(symm_settings.n_atoms * 3 + 1, 0);
    symm_settings.pso_set.max_bound = vector<double>(symm_settings.n_atoms * 3 + 1, 1.0);

    // update function
    auto update_fn = [&](minimize::particle_t &particle) {
        auto sol = particle.get_position();

        sys.set_structure(pso_adapter(sol));
        auto full_structure = sys.get_full_structure();

        sys.update();
        auto value = sys.get_value();

        if (!isnormal(value))
        {
            value = std::numeric_limits<decltype(value)>::max();
        }
        else if (symm_settings.use_cg)
        {
            try {
                minimize::acgsd([&](const auto &pos) {
                    auto structure = sys.get_structure();
                    structure.positions = sp2::dtov3(pos);
                    sys.set_structure(structure);

                    sys.update();
                    return make_pair(sys.get_value(), sys.get_gradient());
                },
                    sp2::v3tod(sys.get_structure().positions),
                    symm_settings.acgsd_set
                );

            } catch (std::domain_error &ex) {}

            bool normal = true;
            for (auto v : sys.get_structure().positions)
                for (auto d : v)
                    normal = normal && isnormal(d);

            if (normal && isnormal(sys.get_value()))
            {
                value = sys.get_value();
                sol = inv_adapter(sys.get_structure());
                full_structure = sys.get_full_structure();
            }
        }

        value += penalty_func(full_structure);

        particle.set_current({value, sol});
    };

    // pso output
    if (symm_settings.pso_set.intermediate_output && comm.rank() == 0)
        io::clear_file(symm_settings.output_filename);

    symm_settings.pso_set.output_fn = [&](const vector<double> &sol) {
        sys.set_structure(pso_adapter(sol));

        // apply symmetries and wrap positions
        auto temp = sys.get_full_structure();

        auto uc = temp.lattice[0][0];
        for (auto &v : temp.positions)
            v -= uc * floor(v / uc);

        // append the structure to the output file
        io::write_structure(symm_settings.output_filename, temp, true);
    };

////////////////////////////////////////////////////////////////////////////////
// Run the PSO                                                                //
////////////////////////////////////////////////////////////////////////////////

    auto sol = minimize::adaptive_pso(update_fn, symm_settings.pso_set, comm);

    // output the final structure
    if (comm.rank() == 0)
    {
        if (!symm_settings.pso_set.intermediate_output)
            io::clear_file(symm_settings.output_filename);

        symm_settings.pso_set.output_fn(sol);
    }

    return EXIT_SUCCESS;
}
