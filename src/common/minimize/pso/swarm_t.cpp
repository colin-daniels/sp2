#include "swarm_t.hpp"

#include "common/util/templates.hpp"
#include "common/util/blas.hpp"
#include "common/util/mpi.hpp"

#include <algorithm>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>

using namespace std;
using namespace sp2;
using namespace minimize;

namespace mpi = boost::mpi;

minimize::swarm_t::swarm_t(pso_update_fn_t update_fn_in,
    const pso_settings_t &settings_in, mpi::communicator comm_in) :
// primary settings
    n_dim(settings_in.min_bound.size()),
    generation(0),
    // current state variables
    evo_factor(0),
    state(evo_state::none),
    // elitist learning
    learning_rate(settings_in.elr_max),
    // best solution found so far
    global_best{numeric_limits<double>::max(), vector<double>(n_dim, 0)},
    // actual input function and mpi communicator
    update_fn(update_fn_in),
    comm(comm_in),
    rng(new util::rng_t()),
    settings(settings_in)
{
    // initialize random number generator
    // with data from /dev/urandom or similar
    rng->seed_random();
}

std::string minimize::swarm_t::get_state() const
{
    switch (state)
    {
    case evo_state::exploration:  return "exploring";
    case evo_state::exploitation: return "exploiting";
    case evo_state::convergence:  return "converging";
    case evo_state::jumping_out:  return "jumping";
    default:
        return "error";
    }
}

double minimize::swarm_t::get_evo_factor() const {
    return evo_factor;}

void minimize::swarm_t::iterate()
{
    if (generation == 0)
        initialize();

    update_evo_factor();
    update_evo_state();

    update_parameters();
    update_els();

    update_particles();
    update_best();

    generation += 1;
}


void minimize::swarm_t::initialize()
{
    // reset particles
    particles.clear();

    for (int i = 0; i < settings.n_particles; ++i)
    {
        // construct
        particles.emplace_back(update_fn,
            settings, rng);

        // randomly generate positions
        auto pos = sol_position_t();
        for (size_t j = 0; j < n_dim; ++j)
            pos.push_back(rng->rand(settings.min_bound[j],
                settings.max_bound[j]));

        // initialize particle solution value/position
        particles[i].set_current({numeric_limits<double>::max(), pos});
    }

    // update particle positions, velocities, and then values
    update_particles();
    // update global_best
    update_best();
}

void minimize::swarm_t::update_evo_factor()
{
    // The evolutionary factor plays into determining the current evolutionary
    // state by being a measure of how far outside the bulk of the swarm
    // the globally best particle is, and is simply defined as:
    //
    //     f = (d_global - d_min) / (d_max - d_min)
    //
    // where the distance d is the average euclidean distance from one particle
    // to all others in the swarm. We approximate d_min and d_max by just using
    // the particles that are the closest and farthest from the average
    // position of all particles in the swarm.

    // get the sum of all particle positions on this process
    auto local_avg = sol_position_t(n_dim, 0);
    for (auto &p : particles)
        vaxpy(1.0, p.get_position(), local_avg);

    // do another sum between all processes
    auto avg_pos = local_avg;
    mpi::all_reduce(comm, local_avg, avg_pos,
        [](const auto &a, const auto &b) {
            auto c = a;
            for (size_t i = 0; i < b.size(); ++i)
                c[i] += b[i];
            return c;
        });

    // recale to get the average
    auto total_particles = settings.n_particles * comm.size();
    vscal(1.0 / total_particles, avg_pos);

    // record the distance from the average position for each particle
    auto avg_dist = vector<double>();
    for (auto &p : particles)
        avg_dist.push_back(vdist(avg_pos, p.get_position()));

    // reset result variables
    dist_min = 0.0;
    dist_max = 0.0;

    // get the best local solution indices
    auto min_idx = nth_index(avg_dist.begin(), avg_dist.end(), 0),
        max_idx = nth_index(avg_dist.begin(), avg_dist.end(), particles.size() - 1);

    // then for each variable, get the best solution from all processes,
    // sum the distance between it and all the other particles, and record it
    for (auto &elem : {
        make_tuple(min_idx,  avg_dist[min_idx], ref(dist_min)),
        make_tuple(max_idx, -avg_dist[max_idx], ref(dist_max))
    }) {
        // get the position corresponding to the particular index and its value
        auto position = particles[get<0>(elem)].get_position();
        auto value = get<1>(elem);

        // find the process with the best one and broadcast its position
        auto best_rank = mpi::all_reduce(comm, value,
            mpi::min_loc<decltype(value)>());
        mpi::broadcast(comm, position, best_rank);

        // get the distance sum between the selected solution and
        // the local particles
        auto dist_sum = 0.0;
        for (auto &p : particles)
            dist_sum += vdist(position, p.get_position());

        // gather the total distance sum via mpi
        auto total_sum = dist_sum;
        mpi::all_reduce(comm, dist_sum, total_sum, plus<decltype(dist_sum)>());

        // set the output variable
        get<2>(elem) = total_sum;
    }

    // get dist_global
    auto temp_dist = 0.0;
    auto global_pos = get<sol_position_t>(global_best);
    for (auto &p : particles)
        temp_dist += vdist(global_pos, p.get_position());

    // sum from all processors
    mpi::all_reduce(comm, temp_dist, dist_global,
        plus<decltype(dist_global)>());

    // correct if we're obviously wrong
    dist_min = min(dist_min, dist_global);
    dist_max = max(dist_max, dist_global);

    dist_min /= total_particles;
    dist_max /= total_particles;
    dist_global /= total_particles;

    // finally, calculate the evolutionary factor
    if (dist_max == dist_min)
        evo_factor = 0.5;
    else
        evo_factor = (dist_global - dist_min) / (dist_max - dist_min);
}

// get state membership [0, 1] for a given evolutionary factor
double get_membership(double evo_factor,
    evo_state state)
{
    // state membership is defined using overlapping piecewise linear functions.
    // membership between the four states is symmetric across evo_factor = 0.5
    // so two cases in this switch statement just fall through
    switch (state)
    {
    case evo_state::exploitation:
        evo_factor = 1 - evo_factor;

        // fallthrough
    case evo_state::exploration:
        if (evo_factor <= 0.4)
            return 0;
        if (evo_factor <= 0.6)
            return 5 * evo_factor - 2;
        if (evo_factor <= 0.7)
            return 1;
        if (evo_factor <= 0.8)
            return -10 * evo_factor + 8;
        return 0;

    case evo_state::jumping_out:
        evo_factor = 1 - evo_factor;

        // fallthrough
    case evo_state::convergence:
        if (evo_factor <= 0.1)
            return 1;
        if (evo_factor <= 0.3)
            return -5 * evo_factor + 1.5;
        return 0;

    default:
        return 0;
    }
}

void minimize::swarm_t::update_evo_state()
{
    auto last_state = state;

    // first get and record the membership values [0, 1] for
    // the current evolutionary factor
    vector<pair<double, evo_state>> states;
    for (auto possible_state : {
        evo_state::exploration,
        evo_state::exploitation,
        evo_state::convergence,
        evo_state::jumping_out
    }) {
        auto mem = get_membership(evo_factor, possible_state);
        // only actually record nonzero ones
        if (mem > 0)
            states.emplace_back(mem, possible_state);
    }

    if (states.empty() || states.size() > 2)
    {
        cerr << states.size() << '\t' << evo_factor << endl;
        throw logic_error("state must be one of the four evolutionary states");
    }
    else if (states.size() == 1)
    {
        // if we are only considered to be in one state
        // no fuzzy logic is needed
        state = get<evo_state>(states.front());
        return;
    }

    // sort states by membership amount
    sort(states.begin(), states.end());

    // this is the first iteration, so just
    // pick the state with the highest membership value
    if (last_state == evo_state::none)
    {
        state = get<evo_state>(states.back());
        return;
    }

    // if any of the states we are a member of _now_ are the state
    // we were in the last iteration, just stay in that state
    // to prevent flip-flopping
    for (auto s : states)
        if (get<evo_state>(s) == last_state)
            return;

    // prefer states that agree with the standard "progression" of states:
    // s1 -> s2 -> s3 -> s4 -> s1 -> ...
    // where 1 is explore, 2 is exploit, 3 is converge, and 4 is jump
    auto theoretical_next = next(last_state);
    for (auto s : states)
    {
        if (get<evo_state>(s) == theoretical_next)
        {
            state = theoretical_next;
            return;
        }
    }

    // all else fails, pick the state with the highest membership
    state = get<evo_state>(states.back());
}

void minimize::swarm_t::update_parameters()
{
    // update per-particle parameters based on the evolutionary factor and state
    for (auto &p : particles)
        p.update_parameters(evo_factor, state);

    // linear decrease in the elitist learning rate based on generation
    learning_rate = settings.elr_max - (settings.elr_max - settings.elr_min)
                                       * (generation / settings.max_generations);
}

void minimize::swarm_t::update_els()
{
    // copy the global best for modification
    auto new_pos = get<sol_position_t>(global_best);

    // random direction to modify
    auto dim = rng->rand(0, n_dim);

    // the dimension will be modified according to a normal distribution with
    // standard deviation = learning_rate
    auto norm_dist = normal_distribution<double>(0.0, learning_rate);

    double dim_val,
        dim_min = settings.min_bound[dim],
        dim_max = settings.max_bound[dim];

    // keep generating a new position for the given dimension until
    // its within the search boundaries
    do {
        dim_val = new_pos[dim] + (dim_max - dim_min)
                                 * norm_dist(rng->get_gen());
    } while (dim_val < dim_min || dim_val > dim_max);

    // update the new position
    new_pos[dim] = dim_val;

    // evaluate it
    particle_t temp(update_fn, settings, rng);
    temp.set_current({numeric_limits<double>::max(), new_pos});
    temp.update();

    auto new_value = temp.get_value();

    // find the best value
    auto best_rank = all_reduce(comm, new_value,
        mpi::min_loc<decltype(new_value)>());

    // get it
    auto best_value = new_value;
    mpi::broadcast(comm, best_value, best_rank);

    // check against current best
    if (isnormal(best_value) &&
        best_value < get<sol_value_t>(global_best))
    {
        // sync/update if its better
        global_best = {new_value, new_pos};
        mpi::broadcast(comm, global_best, best_rank);

        for (auto &p : particles)
            p.update_global(global_best);
    }
    else if (isnormal(new_value))
    {
        // if its not better, replace the particle with the worst fitness
        auto worst_idx = nth_index(particles.begin(), particles.end(),
            particles.size() - 1);

        particles[worst_idx].set_current({new_value, new_pos});
    }
}

void minimize::swarm_t::update_particles()
{
    for (auto &p : particles)
        p.update();
}

void minimize::swarm_t::update_best()
{
    // find the best local solution and update global_best if its better
    auto best_idx = nth_index(particles.begin(), particles.end(), 0);
    global_best = min(global_best, particles[best_idx].get_current());

    // get best solution from all processes
    auto best_rank = mpi::all_reduce(comm, get<sol_value_t>(global_best),
        mpi::min_loc<sol_value_t>());

    mpi::broadcast(comm, global_best, best_rank);

    // update particles
    for (auto &p : particles)
        p.update_global(global_best);
}

std::tuple<double, double, double> minimize::swarm_t::get_dist_info() const
{
    return make_tuple(dist_min, dist_max, dist_global);
}

std::tuple<double, double, double> minimize::swarm_t::get_param_info() const
{
    double inertia,
        accel_self,
        accel_global;

    tie(inertia, accel_self, accel_global) = accumulate(
        particles.begin(), particles.end(), make_tuple(0.0, 0.0, 0.0),
        [](auto sum, auto elem) {
            auto params = elem.get_param_info();
            return make_tuple(
                get<0>(sum) + get<0>(params),
                get<1>(sum) + get<1>(params),
                get<2>(sum) + get<2>(params)
            );
        });

    mpi::all_reduce(comm, double(inertia), inertia, plus<double>());
    mpi::all_reduce(comm, double(accel_self), accel_self, plus<double>());
    mpi::all_reduce(comm, double(accel_global), accel_global, plus<double>());

    auto total_particles = settings.n_particles * comm.size();
    return make_tuple(
        inertia / total_particles,
        accel_self / total_particles,
        accel_global / total_particles
    );
}
