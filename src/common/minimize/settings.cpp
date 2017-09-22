#include "common/minimize/settings.hpp"
#include "common/json/json.hpp"

using namespace std;
using namespace sp2;

bool minimize::ls_settings_t::serialize(Json::Value &output) const
{
    io::serialize_basic(output,
        "armijo_threshold", armijo_threshold,
        "curvature_threshold", curvature_threshold,
        "weak_force_threshold", weak_force_threshold,
        "iteration_limit", iteration_limit,
        "output_level", output_level
    );

    return true;
}

bool minimize::ls_settings_t::deserialize(const Json::Value &input)
{
    io::deserialize_basic(input,
        "armijo_threshold", armijo_threshold,
        "curvature_threshold", curvature_threshold,
        "weak_force_threshold", weak_force_threshold,
        "iteration_limit", iteration_limit,
        "output_level", output_level
    );

    return true;
}

bool minimize::acgsd_settings_t::serialize(Json::Value &output) const
{
    io::serialize_basic(output,
        "linesearch", linesearch_settings,
        "gradient_tolerance", gradient_tolerance,
        "grad_max_tolerance", grad_max_tolerance,
        "value_tolerance", value_tolerance,
        "iteration_limit", iteration_limit,
        "output_level", output_level,
        "except_on_fail", except_on_fail,
        "target_value", target_value,
        "target_ratio_tol", target_ratio_tol,
        "target_exit_min", target_exit_min,
        "intermediate_output_interval", intermediate_output_interval
    );

    return true;
}

bool minimize::acgsd_settings_t::deserialize(const Json::Value &input)
{
    io::deserialize_basic(input,
        "linesearch", linesearch_settings,
        "gradient_tolerance", gradient_tolerance,
        "grad_max_tolerance", grad_max_tolerance,
        "value_tolerance", value_tolerance,
        "iteration_limit", iteration_limit,
        "output_level", output_level,
        "except_on_fail", except_on_fail,
        "target_value", target_value,
        "target_ratio_tol", target_ratio_tol,
        "target_exit_min", target_exit_min,
        "intermediate_output_interval", intermediate_output_interval
    );

    return true;
}

bool minimize::pso_settings_t::serialize(Json::Value &output) const
{
    io::serialize_basic(output,
        "n_particles", n_particles,
        "max_generations", max_generations,
        "intermediate_output", intermediate_output,
        "quiet", quiet);

    io::serialize_basic(output,
        "vel_max_scale", vel_max_scale,
        "accel_delta_min", accel_delta_min,
        "accel_delta_max", accel_delta_max,
        "accel_delta_slight", accel_delta_slight,
        "accel_mag_min", accel_mag_min,
        "accel_mag_max", accel_mag_max,
        "accel_sum_max", accel_sum_max,
        "initial_inertia", initial_inertia,
        "initial_accel_self", initial_accel_self,
        "initial_accel_global", initial_accel_global,
        "elr_max", elr_max,
        "elr_min", elr_min,
        "min_bound", min_bound,
        "max_bound", max_bound
    );

    return true;
}

bool minimize::pso_settings_t::deserialize(const Json::Value &input)
{
    io::deserialize_basic(input,
        "n_particles", n_particles,
        "max_generations", max_generations,
        "intermediate_output", intermediate_output,
        "quiet", quiet);

    io::deserialize_basic(input,
        "vel_max_scale", vel_max_scale,
        "accel_delta_min", accel_delta_min,
        "accel_delta_max", accel_delta_max,
        "accel_delta_slight", accel_delta_slight,
        "accel_mag_min", accel_mag_min,
        "accel_mag_max", accel_mag_max,
        "accel_sum_max", accel_sum_max,
        "initial_inertia", initial_inertia,
        "initial_accel_self", initial_accel_self,
        "initial_accel_global", initial_accel_global,
        "elr_max", elr_max,
        "elr_min", elr_min,
        "min_bound", min_bound,
        "max_bound", max_bound
    );

    return true;
}

bool minimize::fire_settings_t::serialize(Json::Value &output) const
{
    io::serialize_basic(output,
        "grad_tol", grad_tol,
        "max_iter", max_iter,
        "dt_initial", dt_initial,
        "dt_max", dt_max,
        "N_min", N_min,
        "dt_mult_increase", dt_mult_increase,
        "dt_mult_decrease", dt_mult_decrease,
        "alpha_initial", alpha_initial,
        "alpha_mult", alpha_mult
    );

    return true;
}

bool minimize::fire_settings_t::deserialize(const Json::Value &input)
{
    io::deserialize_basic(input,
        "grad_tol", grad_tol,
        "max_iter", max_iter,
        "dt_initial", dt_initial,
        "dt_max", dt_max,
        "N_min", N_min,
        "dt_mult_increase", dt_mult_increase,
        "dt_mult_decrease", dt_mult_decrease,
        "alpha_initial", alpha_initial,
        "alpha_mult", alpha_mult
    );

    return true;
}


bool minimize::metropolis_settings_t::serialize(Json::Value &output) const
{
    io::serialize_basic(output,
        "iteration_limit", iteration_limit,
        "improve_iteration_limit", improve_iteration_limit,
        "output_level", output_level
    );

    return true;
}

bool minimize::metropolis_settings_t::deserialize(const Json::Value &input)
{
    io::deserialize_basic(input,
        "iteration_limit", iteration_limit,
        "improve_iteration_limit", improve_iteration_limit,
        "output_level", output_level
    );

    return true;
}


bool minimize::metropolis_scaling_settings_t::serialize(Json::Value &output) const
{
    io::serialize_basic(output,
        "upscale_by", upscale_by,
        "downscale_by", downscale_by,
        "downscale_max_attempts", downscale_max_attempts
    );

    return true;
}

bool minimize::metropolis_scaling_settings_t::deserialize(const Json::Value &input)
{
    io::deserialize_basic(input,
        "upscale_by", upscale_by,
        "downscale_by", downscale_by,
        "downscale_max_attempts", downscale_max_attempts
    );

    return true;
}

bool minimize::structural_metropolis_settings_t::serialize(Json::Value &output) const
{
    io::serialize_basic(output,
        "enabled", enabled,
        "objective", objective,
        "python_sys_path", python_sys_path,
        "python_module", python_module,
        "python_functions", python_functions,
        "settings", settings
    );

    return true;
}

bool minimize::structural_metropolis_settings_t::deserialize(const Json::Value &input)
{
    io::deserialize_basic(input,
        "enabled", enabled,
        "objective", objective,
        "python_sys_path", python_sys_path,
        "python_module", python_module,
        "python_functions", python_functions,
        "settings", settings
    );

    return true;
}


bool minimize::structural_metropolis_funcs_t::serialize(Json::Value &output) const
{
    io::serialize_basic(output,
        "class", class_,
        "generate", generate,
        "apply", apply,
        "applied", applied,
        "visit", visit,
        "is_repeatable", is_repeatable,
        "scale", scale
    );

    return true;
}

bool minimize::structural_metropolis_funcs_t::deserialize(const Json::Value &input)
{
    io::deserialize_basic(input,
        "class", class_,
        "generate", generate,
        "apply", apply,
        "applied", applied,
        "visit", visit,
        "is_repeatable", is_repeatable,
        "scale", scale
    );

    return true;
}
