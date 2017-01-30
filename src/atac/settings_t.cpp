#include "atac/settings_t.hpp"
#include "common/json/json.hpp"

using namespace sp2;

atac::atac_settings_t::atac_settings_t() :
    iteration_limit(0),
    temperature(1000),
    rate_stw(0),
    rate_divacancy(0),
    rate_add_dimer(0),
    rate_cyclo(0),
    bond_min(0.5),
    bond_max(1.7),
    penalty_scale(500),
    penalty_cutoff(1.7),
    stretch_interval(0),
    stretch_amount{}
{
    min_set.gradient_tolerance = 0;
    min_set.value_tolerance = 1e-4;
    min_set.iteration_limit = 1000;
    min_set.output_level = 0;
    min_set.target_ratio_tol = 5e-3;
    min_set.target_exit_min = 5;
}

bool atac::atac_settings_t::serialize(Json::Value& output) const
{
    io::serialize_basic(output,
        "iteration_limit", iteration_limit,
        "temperature", temperature,
        "rate_stw", rate_stw,
        "rate_divacancy", rate_divacancy,
        "rate_add_dimer", rate_add_dimer,
        "rate_cyclo", rate_cyclo,
        "bond_max", bond_max,
        "bond_min", bond_min,
        "penalty_scale", penalty_scale,
        "penalty_cutoff", penalty_cutoff,
        "stretch_interval", stretch_interval,
        "stretch_amount", stretch_amount,
        "minimize", min_set
    );

    return true;
}

bool atac::atac_settings_t::deserialize(const Json::Value& input)
{
    return io::deserialize_basic(input,
        "iteration_limit", iteration_limit,
        "temperature", temperature,
        "rate_stw", rate_stw,
        "rate_divacancy", rate_divacancy,
        "rate_add_dimer", rate_add_dimer,
        "rate_cyclo", rate_cyclo,
        "bond_max", bond_max,
        "bond_min", bond_min,
        "penalty_scale", penalty_scale,
        "penalty_cutoff", penalty_cutoff,
        "stretch_interval", stretch_interval,
        "stretch_amount", stretch_amount,
        "minimize", min_set
    );
}
