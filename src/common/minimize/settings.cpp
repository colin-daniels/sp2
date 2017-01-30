#include "common/minimize/settings.hpp"
#include "common/json/json.hpp"

using namespace std;
using namespace sp2;

bool minimize::acgsd_settings_t::serialize(Json::Value &output) const
{
    io::serialize_basic(output,
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
