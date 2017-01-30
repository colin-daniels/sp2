#include "particle_t.hpp"
#include "common/util/blas.hpp"

#include <stdexcept>

using namespace std;
using namespace sp2;
using namespace minimize;

minimize::particle_t::particle_t(pso_update_fn_t update_fn_in,
    const pso_settings_t &settings_in,
    std::shared_ptr<util::rng_t> rng_in) :
    inertia(settings_in.initial_inertia),
    accel_self(settings_in.initial_accel_self),
    accel_global(settings_in.initial_accel_global),
    update_fn(update_fn_in),
    rng(rng_in),
    settings(settings_in)
{
    size_t n_dim = settings.min_bound.size();
    auto zero_vec = vector<double>(n_dim, 0);

    current = {numeric_limits<double>::max(), zero_vec};
    self_best   = current;
    global_best = current;

    velocity = zero_vec;
    vel_max  = zero_vec;

    // set velocity maximums to 20% of the search space
    for (size_t i = 0; i < n_dim; ++i)
        vel_max[i] = settings.vel_max_scale
                     * (settings.max_bound[i] - settings.min_bound[i]);
}

const solution_t& minimize::particle_t::get_current() const {
    return current;}

void minimize::particle_t::set_current(const solution_t &sol) {
    current = sol;}

bool minimize::particle_t::operator<(const particle_t &other) const {
    return current < other.current;}


void minimize::particle_t::update_parameters(double evo_factor,
    evo_state state)
{
    update_inertia(evo_factor);
    update_accel(state);
}

void minimize::particle_t::update_global(const solution_t &input) {
    global_best = min(global_best, input);}



void minimize::particle_t::update_inertia(double evo_factor)
{
    // sigmoid mapping for particle inertia
    inertia = 1 / (1 + 1.5 * static_cast<double>(exp(-2.6 * evo_factor)));
}


void minimize::particle_t::update_accel(evo_state state)
{
    // randomly generate delta values
    double delta_ac_self = rng->rand(settings.accel_delta_min,
        settings.accel_delta_max),
        delta_ac_global = rng->rand(settings.accel_delta_min,
        settings.accel_delta_max);

    // determine modification depending on the current state
    switch (state)
    {
    case evo_state::exploration:
        // increase, decrease
        accel_self   += delta_ac_self;
        accel_global -= delta_ac_global;
        break;

    case evo_state::exploitation:
        // slight increase, slight decrease
        accel_self   += settings.accel_delta_slight * delta_ac_self;
        accel_global -= settings.accel_delta_slight * delta_ac_global;
        break;

    case evo_state::convergence:
        // slight increase, slight increase
        accel_self   += settings.accel_delta_slight * delta_ac_self;
        accel_global += settings.accel_delta_slight * delta_ac_global;
        break;

    case evo_state::jumping_out:
        // decrease, increase
        accel_self   -= delta_ac_self;
        accel_global += delta_ac_global;
        break;

    default:
        throw invalid_argument("invalid state");
    }

    // clamp accel values to [1.5, 2.5]
    accel_self   = max(accel_self, settings.accel_mag_min);
    accel_self   = min(accel_self, settings.accel_mag_max);
    accel_global = max(accel_global, settings.accel_mag_min);
    accel_global = min(accel_global, settings.accel_mag_max);

    // normalize values if sum is too large
    if (accel_self + accel_global > settings.accel_sum_max)
    {
        double factor = settings.accel_sum_max / (accel_self + accel_global);
        accel_self   *= factor;
        accel_global *= factor;
    }
}

void minimize::particle_t::update()
{
    auto& position = get<sol_position_t>(current);

    // v_{i + 1} = inertia * v_i
    //     + c_1 rand(0,1) (p_best_i - x_i)
    //     + c_2 rand(0,1) (n_best_i - x_i)
    for (size_t i = 0; i < position.size(); ++i)
    {
        velocity[i] = inertia * velocity[i]
                      + accel_self * rng->rand(0.0, 1.0)
                        * (get<sol_position_t>(self_best)[i] - position[i])
                      + accel_global * rng->rand(0.0, 1.0)
                        * (get<sol_position_t>(global_best)[i] - position[i]);

        // avoid runaway velocity
        velocity[i] = min(velocity[i],  vel_max[i]);
        velocity[i] = max(velocity[i], -vel_max[i]);
        position[i] += velocity[i];

        // check if we're in-bounds (we only update value if we are)
        position[i] = min(position[i], settings.max_bound[i]);
        position[i] = max(position[i], settings.min_bound[i]);
    }

    // reset value and backup
    get<sol_value_t>(current) = numeric_limits<double>::max();
    auto old_sol = current;

    // get the value for the new position
    update_fn(*this);

    if (!std::isnormal(get<sol_value_t>(current)))
        get<sol_value_t>(current) = numeric_limits<double>::max();

    // update historical bests
    self_best = min(self_best, current);
    global_best = min(global_best, current);
}
