#include "common/minimize/minimize.hpp"
#include "swarm_t.hpp"

#include <stdexcept>
#include <iomanip>

using namespace std;
using namespace sp2;
namespace mpi = boost::mpi;

std::vector<double> minimize::adaptive_pso(scalar_fn_t objective_fn,
    const pso_settings_t &settings, boost::mpi::communicator comm)
{
    return adaptive_pso([&](particle_t &particle) {
        auto new_val = objective_fn(particle.get_position());
        particle.set_current({new_val, particle.get_position()});
    }, settings, comm);
}

std::vector<double> minimize::adaptive_pso(pso_update_fn_t update_fn,
    const pso_settings_t &settings, boost::mpi::communicator comm)
{
    if (settings.min_bound.empty() ||
        settings.max_bound.empty() ||
        settings.min_bound.size() != settings.max_bound.size())
        throw std::invalid_argument("Minimum and maximum bounds "
            "must be specified for all dimensions when using the "
            "particle swarm optimization algorithm.");

    // create the swarm object and just iterate until the generation limit is reached
    swarm_t swarm(update_fn, settings, comm);

    auto best = swarm.get_best();
    for (auto gen = 0; gen < settings.max_generations; ++gen)
    {
        swarm.iterate();

        // get the new global best solution (if there is one)
        auto current = swarm.get_best();

        // if we just found a better solution, call the output function
        if (comm.rank() == 0 && settings.intermediate_output && current < best)
            settings.output_fn(get<vector<double>>(current));

        // update the running best solution
        best = min(best, current);

        // output stats for this generation
        if (!settings.quiet)
        {
            int gen_width = ceil(log10(settings.max_generations) + 1);
            double dist_min, dist_max, dist_global,
                inertia, accel_self, accel_global;

            tie(dist_min, dist_max, dist_global) = swarm.get_dist_info();
            tie(inertia, accel_self, accel_global) = swarm.get_param_info();

            if (comm.rank() == 0)
                cout << "g: " << setw(gen_width) << gen << ' '
                     << "v: " << setw(18) << setprecision(14)
                              << get<sol_value_t>(best) << ' '
                     << "f: " << setw(10) << setprecision(6)
                              << swarm.get_evo_factor() << ' '
                     << "s: " << setw(10) << swarm.get_state() << ' '
                     << "d: " << "[" << setw(6)  << setprecision(3)
                              << dist_min << ", "
                              << setw(6)  << setprecision(3)
                              << dist_max << ", "
                              << setw(6)  << setprecision(3)
                              << dist_global << "] "
                     << "a: " << "[" << setw(6)  << setprecision(3)
                              << inertia << ", "
                              << setw(6)  << setprecision(3)
                              << accel_self << ", "
                              << setw(6)  << setprecision(3)
                              << accel_global << "]"
                              << endl;
        }
    }

    return get<sol_position_t>(best);
}
