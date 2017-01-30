#include "atac/system_control_t.hpp"

#include "common/graph/graph.hpp"
#include "common/io/util.hpp"
#include "common/io/file_types/xyz.hpp"

#include "common/minimize/minimize.hpp"
#include "common/structure_t.hpp"

#include <unordered_map>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>
#include <tuple>

#include <mpi.h>

using namespace std;
using namespace sp2;

/// mutation error type
const atac::mut_instance_t mut_error{"error", 0, [](){return false;}, false};

atac::system_control_t::system_control_t(util::mpi_group_t mpi_group) :
    iteration(0), accepted_mutations(0), best_delta(numeric_limits<double>::max()),
    mpi_info(mpi_group), sub_sys(nullptr)
{
    int seed = rand();
    MPI_Bcast(&seed, 1, MPI_INT, 0, mpi_info.group_comm);
    srand(seed);
}

void atac::system_control_t::bcast(MPI_Comm comm, int root)
{
    // get a temporary of the current structure
    structure_t structure_copy = get_structure();
    // sync it with all other MPI tasks
    structure_copy.bcast(comm, root);
    // set our own structure to the sync'd value
    set_structure(structure_copy);

    // send bond graph and potentials
    graph.bcast(comm, root);

    MPI_Bcast(&total_potential, 1, MPI_DOUBLE, root, comm);
    MPI_Bcast(&penalty_potential, 1, MPI_DOUBLE, root, comm);
}

void atac::system_control_t::set_structure(const structure_t &input)
{
    copy_n(input.lattice[0], 9, lattice[0]);
    type = input.types;
    position = input.positions;
    gradient.resize(position.size(), 0);
    bond_control.set_lattice(lattice);

    na = position.size() / 3;
}

structure_t atac::system_control_t::get_structure() const
{
    structure_t output;

    copy_n(lattice[0], 9, output.lattice[0]);
    output.types = type;
    output.positions = position;

    return output;
}

bool atac::system_control_t::stretch_lattice()
{
    // modify system by making the lattice vectors longer (or shorter)
    for (int i = 0; i < 3; ++i)
    {
        double amount = settings.stretch_amount[i];
        if (amount == 0)
            continue;

        double len = 0;
        for (int j = 0; j < 3; ++j)
            len += lattice[i][j] * lattice[i][j];

        len = sqrt(len);

        // modify in the direction of the actual lattice vector
        // (not necessarily just x, y, or z)
        for (int j = 0; j < 3; ++j)
            lattice[i][j] += amount * lattice[i][j] / len;
    }

    return true;
}


void atac::system_control_t::init(std::unique_ptr<atac_subsys_t> &&sys_in,
    atac_settings_t input_settings)
{
    // copy settings
    settings = input_settings;

    // reset iterations
    iteration = mpi_info.group_num;
    accepted_mutations = 0;

    // register mutations by adding structures containing
    // their names, probabilities, and application functions
    // (note the value capture of the 'this' pointer in the lambda functions)
    mutations = {
        // Stone-Thrower-Wales rotation ////////////////////////////////
        {"stw", settings.rate_stw,       [=](){return stone_wales();},   false},
        // two atom vacancy (divacancy) ////////////////////////////////
        {"div", settings.rate_divacancy, [=](){return divacancy();},     false},
        // two atom add dimer //////////////////////////////////////////
        {"add", settings.rate_add_dimer, [=](){return add_dimer();},     false},
        // cycloaddition (sp2->sp3->sp2) ///////////////////////////////
        {"cyc", settings.rate_cyclo,     [=](){return cycloaddition();},  true}
    };

    // take ownership of the subsystem object passed in
    sub_sys = move(sys_in);
    auto initial_structure = sub_sys->get_structure();

    // setup the bond controller
    bond_control.init(initial_structure.lattice, settings.bond_max, 0);

    // load position, type, and lattice data
    set_structure(initial_structure);
    na = position.size() / 3;

    bond_control.update(position);

    // and copy bonds which are greater than bond_min
    graph::ud_graph_t bc_graph = bond_control.get_graph();
    vector<double> length = bond_control.get_bond_lengths();

    vector<pair<int, int>> bond_pairs;
    for (auto edge : bc_graph.edges())
        if (length[edge.id] > settings.bond_min)
            bond_pairs.push_back(edge);

    // initialize the actual bond graph with the 'good' bonds
    graph = graph::ud_graph_t(bond_pairs, na);

    // do initial force/potential calculation
    update();
}

void atac::system_control_t::update()
{
    // update simple counters
    na = type.size();
    nb = graph.n_edges();

    // set subsystem structural data
    sub_sys->set_structure(get_structure());

    // perform the system update and copy values back
    sub_sys->update();

    total_potential = sub_sys->get_value();
    gradient = sub_sys->get_gradient();

    // calculate the penalty potential
    update_penalty();
}

void atac::system_control_t::update_penalty()
{
    // reset
    penalty_potential = 0;

    // get the (x,y,z) components of the bonds defined in the system graph
    auto deltas = bond_control.calc_delta_graph(graph, position);

    // for all enforced bonds
    for (auto edge : graph.edges())
    {
        double dist = 0;
        for (int k = 0; k < 3; ++k)
            dist += deltas[edge.id * 3 + k] * deltas[edge.id * 3 + k];

        // only apply the penalty potential if atoms that
        // are forced to be adjacent are farther away
        // than the cutoff distance
        if (dist <= settings.penalty_cutoff * settings.penalty_cutoff)
            continue;

        dist = sqrt(dist);

        // the potential itself is just a quadratic
        double d_cut = dist - settings.penalty_cutoff;
        penalty_potential += settings.penalty_scale * d_cut * d_cut;

        double grad_mag = 2 * settings.penalty_scale * d_cut / dist;
        for (int k = 0; k < 3; ++k)
        {
            gradient[edge.a * 3 + k] += grad_mag * -deltas[edge.id * 3 + k];
            gradient[edge.b * 3 + k] -= grad_mag * -deltas[edge.id * 3 + k];
        }
    }

    // update the total
    total_potential += penalty_potential;
}

bool atac::system_control_t::iterate()
{
    if (iteration >= settings.iteration_limit)
    {
        if (mpi_info.world_rank == 0)
            cout << "ATAC Iteration limit reached. "
                 << accepted_mutations << " of " << iteration
                 << " mutations accepted." << endl;
        return false;
    }

    // backup
    double backup_potential = total_potential,
        backup_penalty = penalty_potential;
    graph::ud_graph_t backup_graph = graph;
    structure_t backup_structure = get_structure();

    // mutate
    mut_instance_t mut;

    if (settings.stretch_interval > 0 &&
        (iteration + 1) % settings.stretch_interval == 0)
    {
        mut.name = "lat";
        mut.always_accept = true;
        stretch_lattice();
    }
    else
    {
        // randomly pick and enact a mutation from our mutation set
        mut = attempt_mutation();
    }

    if (mut.name == "error")
    {
        cout << "Error, mutation attempts exceeded maximum allowed."
             << endl;
        return false;
    }

    // calculate target energy, aka the energy which the system must be below
    // for it to be accepted by Arrhenius' equation (P = e^(-delta E / kT))
    // (note: just generate P and solve for delta E)
    double target = total_potential - (boltzmann * settings.temperature)
                                      * log(util::rand_double(0, 1));

    // minimization, updates total_potential and relaxes the structure
    if (mut.always_accept)
    {
        // disable targeted minimization if we are going to always
        // accept this mutation
        auto modified_settings = settings.min_set;

        modified_settings.target_ratio_tol = 0;
        target = total_potential;

        minimize::acgsd(*this, modified_settings);
    }
    else
    {
        settings.min_set.target_value = target;
        minimize::acgsd(*this, settings.min_set);
    }

    // a mutation is accepted if the energy is not infinite or NaN,
    // the new energy is less than the target energy, and the physical potential
    // is also less than the physical potential part of the target energy
    int mut_accepted = isfinite(total_potential) &&
                       (mut.always_accept  ||
                        (total_potential < target &&
                         total_potential - penalty_potential < target - backup_penalty));

    vector<int> accepted_list(mpi_info.world_size, 0);
    MPI_Allgather(&mut_accepted, 1, MPI_INT,
        accepted_list.data(), 1, MPI_INT, mpi_info.world_comm);

    auto accepted_loc = find(accepted_list.begin(),
        accepted_list.end(), 1);

    // if nobody accepted a mutation
    // (aka, no elemented of accepted_list were nonzero)
    if (accepted_loc == accepted_list.end())
    {
        // update the best delta since the last accepted iteration
        double delta = total_potential - target;
        if (delta > best_delta || !isfinite(delta))
            delta = best_delta;

        MPI_Allreduce(&delta, &best_delta, 1, MPI_DOUBLE, MPI_MIN,
            mpi_info.world_comm);

        if (mpi_info.world_rank == 0)
            cout << setw(6) << accepted_mutations << " "
                 << setw(6) << iteration << " "
                 << "reject "
                 << "target_delta: " << setprecision(6) << setw(12) << best_delta
                 << endl;

        // reset
        total_potential = backup_potential;
        penalty_potential = backup_penalty;
        graph = backup_graph;
        set_structure(backup_structure);

        iteration += mpi_info.num_groups;
    }
    else
    {
        int accepted_rank = distance(accepted_list.begin(), accepted_loc);

        // update others with the new structure
        bcast(mpi_info.world_comm, accepted_rank);

        // output
        if (mpi_info.world_rank == accepted_rank)
        {
            cout << setw(6) << accepted_mutations << " "
                 << setw(6) << iteration << " "
                 << "accept "
                 << mut.name << " "
                 << "energy: " << setprecision(6) << setw(12)
                 << total_potential - penalty_potential << " "
                 << "delta: " << setprecision(6) << setw(12)
                 << total_potential - backup_potential << endl;

            sub_sys->append_output("mutate.xyz");
        }

        // increment number of accepted mutations
        accepted_mutations += 1;

        // increment number of iterations up to the one that was found
        int iter_increment = mpi_info.group_num + 1;
        MPI_Bcast(&iter_increment, 1, MPI_INT, accepted_rank,
            mpi_info.world_comm);

        iteration += iter_increment;

        // reset best delta when something is accepted (since the structure changes)
        best_delta = numeric_limits<double>::max();
    }

    return true;
}

atac::mut_instance_t atac::system_control_t::attempt_mutation(int max_attempts)
{
    for (int attempts = 0; attempts < max_attempts; ++attempts)
    {
        mut_instance_t mut = pick_mutation();

        // if it applies successfully, return it
        if (mut.apply())
            return mut;

        cout << mut.name << endl;
    }

    // if we exceed the maximum number of mutation attempts, return an error
    return mut_error;
}

atac::mut_instance_t atac::system_control_t::pick_mutation()
{
    // generate a random number to pick the mutation
    double rand_n = util::rand_double(0, 1);

    // go through the mutations
    for (auto mut : mutations)
    {
        // if the random # is less than our probability, return the type
        if (rand_n < mut.probability)
            return mut;

        // otherwise, subtract off the probability 'so far' and keep going
        rand_n -= mut.probability;
    }

    // all probabilities should add up to one, so we shouldn't get here normally
    return mut_error;
}

std::vector<graph::ud_edge_t>
atac::system_control_t::random_sp2_carbon_bonds()
{
    vector<graph::ud_edge_t> edges;
    for (auto edge : graph.edges())
        if (graph.degree(edge.a) == 3 &&
            graph.degree(edge.b) == 3 &&
            type[edge.a] == atom_type::CARBON &&
            type[edge.b] == atom_type::CARBON)
            edges.push_back(edge);

    // randomize order
    random_shuffle(edges.begin(), edges.end());
    return edges;
}

bool atac::system_control_t::stone_wales()
{
    int err_size = 0,
        err_adj  = 0;

    for (auto edge : random_sp2_carbon_bonds())
    {
        int a = edge.a,
            b = edge.b;

        vector<int> neigh_a, neigh_b;
        for (int id : graph.neighbors(a))
            if (id != b) neigh_a.push_back(id);

        for (int id : graph.neighbors(b))
            if (id != a) neigh_b.push_back(id);

        if ((neigh_a.size() != 2 || neigh_b.size() != 2) && ++err_size)
            continue;

        // check for shared (or close) neighbors
        // to avoid making triangles and squares
        auto are_close = [&](int c, int d) {
            if (c == d || graph.adjacent(c, d))
                return true;

            for (auto id : graph.neighbors(c))
                if (id != a && id != b && graph.adjacent(d, id))
                    return true;

            return false;
        };

        if ((are_close(neigh_a[0], neigh_b[0]) ||
             are_close(neigh_a[0], neigh_b[1]) ||
             are_close(neigh_a[1], neigh_b[1])) && ++err_adj)
            continue;

        // need to swap randomly since neighbors are probably sorted
        if (rand() % 2)
            swap(neigh_a[0], neigh_a[1]);

        // configure the ids so that neigh_a[0] and neigh_b[0]
        // are the farthest away from one another
        auto get_dist = [&](int c, int d) {
            return bond_control.calc_distance(&position[c * 3], &position[d * 3]);
        };

        double dist_a0b0 = get_dist(neigh_a[0], neigh_b[0]),
            dist_a0b1 = get_dist(neigh_a[0], neigh_b[1]);

        if (dist_a0b0 < dist_a0b1)
            swap(neigh_b[0], neigh_b[1]);

        // modify the graph (do the rotation)
        graph.remove_edge(a, neigh_a[0]);
        graph.remove_edge(b, neigh_b[0]);

        graph.add_edge(a, neigh_b[0]);
        graph.add_edge(b, neigh_a[0]);

        // reposition atoms to an initial guess (could be done better)
        // of what the structure will be like after relaxation
        double delta_a[3] = {};
        bond_control.calc_delta(&position[neigh_a[1] * 3],
            &position[neigh_b[0] * 3], delta_a);

        double delta_b[3] = {};
        bond_control.calc_delta(&position[neigh_a[0] * 3],
            &position[neigh_b[1] * 3], delta_b);

        for (int k = 0; k < 3; ++k)
        {
            position[a * 3 + k] = position[neigh_a[1] * 3 + k] + 0.5 * delta_a[k];
            position[b * 3 + k] = position[neigh_a[0] * 3 + k] + 0.5 * delta_b[k];
        }

        // done
        return true;
    }

    cout << "stw failed: " << err_adj << '\t' << err_size << endl;

    return false;
}

bool atac::system_control_t::divacancy()
{
    for (auto edge : random_sp2_carbon_bonds())
    {
        int a = edge.a,
            b = edge.b;

        vector<int> neigh_a, neigh_b;
        for (int id : graph.neighbors(a))
            if (id != b) neigh_a.push_back(id);

        for (int id : graph.neighbors(b))
            if (id != a) neigh_b.push_back(id);

        // check for situations that would make duplicate bonds
        if (graph.adjacent(neigh_a[0], neigh_a[1]) ||
            graph.adjacent(neigh_b[0], neigh_b[1]))
            continue;

        graph.add_edge(neigh_a[0], neigh_a[1]);
        graph.add_edge(neigh_b[0], neigh_b[1]);

        // removal invalidates vertex ids larger than the input,
        // so make sure we remove in the right order
        if (a < b) swap(a, b);

        remove_atom(a);
        remove_atom(b);

        // done
        return true;
    }

    return false;
}

bool atac::system_control_t::add_dimer()
{
    // TODO: add_dimer
    return true;
}

bool atac::system_control_t::cycloaddition()
{
    auto edges = random_sp2_carbon_bonds();

    for (auto edge_1 : edges)
    {
        int a = edge_1.a,
            b = edge_1.b;

        for (auto edge_2 : edges)
        {
            int c = edge_2.a,
                d = edge_2.b;

            if (edge_1.id >= edge_2.id)
                continue;

            // basically just make sure the atoms aren't near eachother
            // in terms of bonds (not real space)
            int g_dist = graph::bfs_dist(graph, a, c, 12);
            if (g_dist != -1)
                continue;

            auto get_dist = [&](double *pos_a, double *pos_b) {
                return bond_control.calc_distance(pos_a, pos_b);};

            double dist_ac = get_dist(&position[a * 3], &position[c * 3]),
                dist_ad = get_dist(&position[a * 3], &position[d * 3]),
                dist_bc = get_dist(&position[b * 3], &position[c * 3]),
                dist_bd = get_dist(&position[b * 3], &position[d * 3]);

            if ((dist_ac + dist_bd) > (dist_bc + dist_ad))
                swap(a, b);

            // distance cutoff
            if (max(dist_ac, dist_bd) > 2.5)
                continue;

            // sp2 -> sp3
            graph.add_edge(a, c);
            graph.add_edge(b, d);

            // and then back from sp3 -> sp2
            graph.remove_edge(a, b);
            graph.remove_edge(c, d);

            // done
            return true;
        }
    }

    return false;
}

void atac::system_control_t::remove_atom(int id)
{
    // get system size to be safe
    na = type.size();
    if (id >= na)
        return; // out of bounds

    // swap all the data with the last atom (to avoid
    // having to move all the others down one)
    swap(type[id], type[na - 1]);
    for (int i = 0; i < 3; ++i)
    {
        swap(position[id * 3 + i], position[(na - 1) * 3 + i]);
        swap(gradient[id * 3 + i], gradient[(na - 1) * 3 + i]);
    }

    // decrement the atom counter
    na -= 1;

    type.resize(na);
    position.resize(na * 3);
    gradient.resize(na * 3);

    graph.remove_vertex(id);
}

int atac::system_control_t::create_atom(atom_type type_in, double *pos)
{
    na = type.size() + 1;

    type.push_back(type_in);
    for (int i = 0; i < 3; ++i)
    {
        position.push_back(pos[i]);
        gradient.push_back(0.0);
    }

    return graph.add_vertex();
}
