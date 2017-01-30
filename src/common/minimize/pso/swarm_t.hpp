#ifndef SP2_SWARM_T_HPP
#define SP2_SWARM_T_HPP

#include "common/util/random.hpp"
#include "common/minimize/minimize.hpp"
#include "particle_t.hpp"

#include <vector>
#include <memory>
#include <tuple>

#include <boost/mpi/communicator.hpp>

namespace sp2 {
namespace minimize {

/// particle swarm class for the adaptive particle swarm optimization algorithm
class swarm_t
{
private:
    /// the dimensionality of the search space
    size_t n_dim,
    /// current generation
        generation;

    /// min average euclidean distance from one particle to all others
    double dist_max,
    /// max average euclidean distance from one particle to all others
        dist_min,
    /// distance from the global best particle to all others
        dist_global,
    /// current evolutionary factor
        evo_factor;
    /// current evolutionary state
    evo_state state;
    /// current elitist learning rate
    double learning_rate;

    /// particles belonging to this process
    std::vector<particle_t> particles;
    /// current best solution
    solution_t global_best;

    /// objective function used to calculate solution value
    pso_update_fn_t update_fn;
    /// mpi communicator that the swarm operates on
    boost::mpi::communicator comm;
    /// rng object used for all random number generation
    std::shared_ptr<util::rng_t> rng;

    /// pso settings
    pso_settings_t settings;
public:
    swarm_t(pso_update_fn_t update_fn_in, const pso_settings_t &settings_in,
        boost::mpi::communicator comm_in);

    /// perform one iteration/simulate one generation of the swarm
    void iterate();
    /// get the current best solution to the objective function
    solution_t get_best() const {return global_best;}
    /// get the current evolutionary state string
    std::string get_state() const;
    /// get the current evolutionary factor
    double get_evo_factor() const;
    /// get dist_min, dist_max, dist_global
    std::tuple<double, double, double> get_dist_info() const;
    /// get average inertia, accel_self, accel_global
    std::tuple<double, double, double> get_param_info() const;

private:
    /// initialize/reset the system
    void initialize();

    /// calculate the evolutionary factor
    void update_evo_factor();
    /// calculate the current evolutionary state
    void update_evo_state();
    /// update particle parameters based on current evolutionary state
    void update_parameters();
    /// perform the elitest learning strategy for global_best
    void update_els();
    /// update positions, velocities, and values for all particles
    void update_particles();
    /// update global best
    void update_best();
};

} // namespace minimize
} // namespace sp2

#endif // SP2_SWARM_T_HPP
