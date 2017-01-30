#ifndef SP2_PARTICLE_T_HPP
#define SP2_PARTICLE_T_HPP

#include "common/util/random.hpp"
#include "common/minimize/minimize.hpp"
#include "common/minimize/settings.hpp"

#include <vector>
#include <memory>
#include <tuple>

namespace sp2 {
namespace minimize {

/// swarm evolutionary state
enum class evo_state : int
{
    none         = -1,
    exploration  =  0,
    exploitation =  1,
    convergence  =  2,
    jumping_out  =  3
};

/// return the "normal" next state for an input state
static inline evo_state next(evo_state state) {
    return static_cast<evo_state>((static_cast<int>(state) + 1) % 4);}


// the "solution types" representing the value and position of a given
using sol_value_t    = double;
using sol_position_t = std::vector<double>;
using solution_t     = std::pair<sol_value_t, sol_position_t>;


/// single particle for particle swarm optimization
class particle_t
{
private:
    /// inertia factor (omega)
    double inertia,
    /// scale factor for acceleration towards self_best
        accel_self,
    /// scale factor for acceleration towards global_best
        accel_global;

    /// best solution found by this particle (historically)
    solution_t self_best,
    /// best solution found by any particle (historically)
        global_best,
    /// current solution
        current;

    /// current particle velocity
    std::vector<double> velocity,
    /// velocity maximum for each dimension
        vel_max;

    /// objective function used to update value for the current position
    pso_update_fn_t update_fn;
    /// rng object used to generate all random numbers
    std::shared_ptr<util::rng_t> rng;

    /// pso settings
    pso_settings_t settings;
public:
    particle_t(pso_update_fn_t update_fn_in,
        const pso_settings_t &settings_in,
        std::shared_ptr<util::rng_t> rng_in = std::make_shared<util::rng_t>());

    /// get the current solution
    const solution_t& get_current() const;
    /// set the current solution
    void set_current(const solution_t &sol);

    /// update parameters (accel, inertia)
    void update_parameters(double evo_factor, evo_state state);
    /// update position, velocity, and current value
    void update();
    /// update global_best
    void update_global(const solution_t &input);

    /// comparison operator for sorting
    bool operator<(const particle_t &other) const;

    /// get current value
    sol_value_t get_value() const {
        return get<sol_value_t>(current);}

    /// get current position
    const sol_position_t& get_position() const {
        return get<sol_position_t>(current);}

    /// get parameters: inertia, accel_self, accel_global
    std::tuple<double, double, double> get_param_info() const {
        return std::make_tuple(inertia, accel_self, accel_global);}
private:
    void update_inertia(double evo_factor);
    void update_accel(evo_state state);
};

} // namespace minimize
} // namespace sp2

#endif // SP2_PARTICLE_T_HPP
