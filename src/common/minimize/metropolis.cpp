#include "common/minimize/minimize.hpp"

#include <vector>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace sp2;

std::vector<double> minimize::metropolis(
    scalar_fn_t objective_fn,
    mutation_fn_t mutation_fn,
    std::vector<double> initial_position,
    const metropolis_settings_t &settings)
{
    if (settings.iteration_limit <= 0 &&
        settings.improve_iteration_limit <= 0)
        throw invalid_argument("all exit conditions disabled or set to 0.");

    int mutations = 0,        // total number of mutations tested
        mutations_since = 0;  // number of mutations since last improvement

    // Initialize
    vector<double> position_best = initial_position;
    double value_best = objective_fn(position_best);

    if (settings.output_level > 1) {
        cout << " Initial Value: " << value_best << endl;
    }

    do {
        auto position_new = position_best;

        mutation_fn(position_new);
        mutations += 1;
        mutations_since += 1;

        double value_new = objective_fn(position_new);

        if (value_new < value_best) {
            if (settings.output_level > 1) {
                cout << "Improved Value: " << value_new
                     << "   after " << mutations_since << " mutations"
                     << endl;
            }

            // keep modified position and value
            value_best = value_new;
            position_best = std::move(position_new);
            mutations_since = 0;
        }

        // for now, simply reject all inferior values
        // (this is Metropolis at 0K)

        // exit conditions
        bool done = false;
        done = done || settings.improve_iteration_limit > 0 &&
                       mutations_since >= settings.improve_iteration_limit;

        done = done || settings.iteration_limit > 0 &&
                       mutations >= settings.iteration_limit;

        if (done)
        {
            // output final info
            if (settings.output_level > 1)
                cout << "Metropolis Finished.\n"
                     << "     Value: " << value_best << '\n'
                     << " Mutations: " << mutations << '\n'
                     << endl;

            break;
        }
    } while (true);

    return position_best;
}
