#ifndef SP2_MINIMIZE_SETTINGS_HPP
#define SP2_MINIMIZE_SETTINGS_HPP

#include <vector>
#include <functional>

#include "common/json/json_serializable_t.hpp"
#include "common/function_types.hpp"

namespace sp2 {
namespace minimize {

////////////////////////////////////////////////////////////////////////////////
// Input settings structures for the minimization routines.                   //
////////////////////////////////////////////////////////////////////////////////

/// output function type, takes a vector and does
/// something with it (presumably outputs to a file)
using output_fn_t = void_fn_t<const std::vector<double> &>;

/// input settings structure for conjugate gradient minimization
struct acgsd_settings_t : public io::json_serializable_t
{
    /// exit condition: tol > ||gradient||, 0 to disable
    double gradient_tolerance = 1e-4;
    /// exit condition: tol > ||gradient||_{max norm}, 0 to disable
    double grad_max_tolerance = 0;
    /// exit condition: tol > |delta value|, 0 to disable
    double value_tolerance = 0;

    /// exit condition: limit > iteration, 0 to disable
    int iteration_limit = 1000;

    /// how much output should be sent to standard output (0-3):
    /// 0: None
    /// 1: Errors only
    /// 2: Final information after minimization
    /// 3: Value and forces every iteration
    int output_level = 3;

    /// whether or not an exception is thrown when linesearch fails
    bool except_on_fail = false;

    /// value for targeted minimization
    double target_value = 0;

    /// exit condition for targeted minimization:
    /// |delta value / (value - target)| < target_ratio_tol,
    /// 0 to disable
    double target_ratio_tol = 0;

    /// number of iterations in a row that must
    /// fail the target ratio exit condition
    /// test before minimization is aborted,
    /// 0 to disable
    int target_exit_min = 0;

    /// how often to call the output function with the current
    /// state, set to 0 for never
    int intermediate_output_interval = 0;

    /// function to call with the current state
    /// every intermediate_output_interval iterations
    output_fn_t output_fn = [](const auto &) {};

    bool serialize(Json::Value &output) const;

    bool deserialize(const Json::Value &input);
};


/// simple settings structure for APSO
struct pso_settings_t : public io::json_serializable_t
{
    /// number of particles per PROCESS
    int n_particles = 100,
    /// maximum number of generations (exit condition)
        max_generations = 1000;

    /// factor used when calculating velocity maxima for each dimension,
    /// i.e. a value of 0.2 will make the maximum velocity equal to 20%
    /// of the search space range for a given dimension
    double vel_max_scale = 0.2,
    /// minimum change in acceleration values for each generation
        accel_delta_min = 0.05,
    /// maximum change in acceleration values for each generation
        accel_delta_max = 0.1,
    /// scaling factor for "slight" accel changes for each generation
        accel_delta_slight = 0.5,
    /// minimum magnitude for each individual accel coefficient
        accel_mag_min = 1.5,
    /// maximum magnitude for each individual accel coefficient
        accel_mag_max = 2.5,
    /// maximum bound on the sum of the two accel coefficients (self + global)
        accel_sum_max = 4.0;

    /// initial inertia factor (omega)
    double initial_inertia = 0.9,
    /// initial scale factor for acceleration towards self_best
        initial_accel_self = 2.0,
    /// initial scale factor for acceleration towards global_best
        initial_accel_global = 2.0;

    // maximum elitist learning rate
    double elr_max = 1.0;
    /// minimum elitist learning rate
    double elr_min = 0.1;

    /// lower bound of search space for each dimension, required
    std::vector<double> min_bound,
    /// upper bound of search space for each dimension, required
        max_bound;

    /// whether to call the output function every
    /// time a better solution is found
    bool intermediate_output = true;

    /// whether to output anything to console
    bool quiet = false;

    /// output function to call when a better solution is found
    /// (assuming intermediate_output is set to true)
    output_fn_t output_fn = [](const auto &) {};

    bool serialize(Json::Value &output) const;

    bool deserialize(const Json::Value &input);
};

} // namespace minimize
} // namespace sp2

#endif // SP2_MINIMIZE_SETTINGS_HPP