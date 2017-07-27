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


struct fire_settings_t : public io::json_serializable_t
{
    /// number of iterations we have been "going downhill"
    double N_min = 5,
    /// amount to increase dt when going downhill
        dt_mult_increase = 1.1,
    /// amount to decrease dt when going uphill
        dt_mult_decrease = 0.5,
    /// opposite of inertia
        alpha_initial = 0.1,
    /// amount to decrease alpha if going downhill
        alpha_mult = 0.99;

    /// initial time step
    double dt_initial = 0.001,
    /// max time step
        dt_max = 0.05;

    /// maximum number of iterations
    int max_iter = 1000;
    /// gradient tolerance exit condition
    double grad_tol = 1e-5;

    /// function to call with the current state
    /// every intermediate_output_interval iterations
    output_fn_t output_fn = [](const auto &) {};

    bool serialize(Json::Value &output) const;

    bool deserialize(const Json::Value &input);
};


/// input settings structure for metropolis minimization
struct metropolis_settings_t : public io::json_serializable_t
{
    /// exit condition: total number of mutations tested, 0 to disable
    int iteration_limit = 0;
    /// exit condition: number of mutations tested since last improvement, 0 to disable
    int improve_iteration_limit = 100;

    /// how much output should be sent to standard output (0-3):
    /// 0: None
    /// 1: Errors only
    /// 2: Final information after minimization
    /// 3: Information after each successful refinement
    int output_level = 3;

    bool serialize(Json::Value &output) const;

    bool deserialize(const Json::Value &input);
};


struct structural_metropolis_funcs_t : public io::json_serializable_t
{
    /// When false, "basic callbacks" are used.
    /// When true, "advanced callbacks" are used.
    bool advanced = false;

    /// Basic callback for an all-in-one mutation function.
    ///
    /// All arguments supplied to this are keyword arguments.
    /// They are not currently documented (sorry).
    /// For future compatibility, the callback should always take a **kw
    ///  argument to collect any remainder not explicitly of interest.
    ///
    /// This should return a value which can be recognized
    ///  as a 'structural_mutation_t'.
    std::string mutate = "mutate";

    /// Advanced callback for generating mutations.
    ///
    /// If defined, then it should accept the same arguments
    /// as 'mutate', but can return any arbitrary python object.
    /// The other callbacks will decide how to interpret this object.
    ///
    /// Documentation for other advanced callbacks will use the term "mutation"
    /// to indicate whatever form of data is returned by this function.
    std::string generate = "generate";

    /// Advanced callback for applying a mutation to a structure.
    ///
    /// If defined, then it should accept:
    ///  (1) a mutation,
    ///  (2) the standard kw arguments provided to 'mutate'
    /// It should produce a value parsable as 'structural_mutation_t'.
    ///
    /// If omitted, will search for a function named 'apply',
    /// and then if that is not found, the following definition is assumed:
    ///
    /// def apply(mutation, **kw):
    ///     return mutation  # assume mutation is 'structural_mutation_t'
    std::string apply;

    /// Advanced callback for signaling successful mutations. (lower energy)
    ///
    /// If defined, then it should accept:
    ///  (1) a mutation,
    ///  (2) the standard kw arguments provided to 'mutate'
    /// It should produce 'None'.
    ///
    /// This function will be called when a mutation is accepted by the
    /// metropolis algorithm.
    ///
    /// If omitted, will search for a function named 'accept',
    /// and then if that is not found, the following definition is assumed:
    ///
    /// def accept(mutation, **kw):
    ///     pass
    std::string accept;

    /// Advanced callback for determining if mutations are worth repeating.
    /// (i.e., if the mutation decreases total energy, might we expect that
    ///        repeating the mutation decreases energy again?)
    ///
    /// If defined, then it should accept:
    ///  (1) a mutation,
    ///  (2) the standard kw arguments provided to 'mutate'
    /// It should produce a 'bool'
    ///
    /// If omitted, will search for a function named 'is_repeatable',
    /// and then if that is not found, the following definition is assumed:
    ///
    /// def is_repeatable(mutation, **kw):
    ///     return false  # never repeat anything.
    std::string is_repeatable;

    /// Advanced callback for modifying a mutation's strength.
    ///
    /// If defined, then it should accept:
    ///  (1) a mutation,
    ///  (2) a float 'strength' factor (< 1 or > 1),
    ///  (3) the standard kw arguments provided to 'mutate'
    /// It should produce a mutation.
    ///
    /// If omitted, will search for a function named 'scale',
    /// and then if that is not found, the following definition is assumed:
    ///
    /// def scale(mutation, factor, **kw):
    ///     return mutation
    std::string scale;

    virtual bool serialize(Json::Value &output) const;
    virtual bool deserialize(const Json::Value &input);
};


struct metropolis_scaling_settings_t : public io::json_serializable_t
{
    double upscale_by = 1.1;
    double downscale_by = 0.5;
    int downscale_max_attempts = 3;

    bool serialize(Json::Value &output) const;

    bool deserialize(const Json::Value &input);
};

struct structural_metropolis_settings_t : public io::json_serializable_t
{
    bool enabled = false;

    // Directories to be prepended to sys.path, where python modules may be found.
    //
    // The default prepends an entry of "" (the current directory), just like
    // the python interpreter itself typically does behind the scenes.
    std::vector<std::string> python_sys_path = {""};
    std::string python_module = "mutate";
    structural_metropolis_funcs_t python_functions;
    minimize::metropolis_settings_t settings;
    minimize::metropolis_scaling_settings_t scaling_settings;

    virtual bool serialize(Json::Value &output) const;
    virtual bool deserialize(const Json::Value &input);
};


} // namespace minimize
} // namespace sp2

#endif // SP2_MINIMIZE_SETTINGS_HPP
