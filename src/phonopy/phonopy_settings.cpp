//
// Created by cc on 11/26/16.
//

#include "phonopy_settings.hpp"
#include "json/json.h"
#include "common/json/json.hpp"
#include <utility>

using namespace std;

std::string check_pol_axes(const Json::Value &input)
{
    if (!input || !input.isString())
        return "";

    return input.asString();
}

bool sp2::phonopy::phonopy_settings_t::serialize(Json::Value &output) const
{
    io::serialize_basic(output,
        "n_samples", n_samples,
        "minimize", min_set,
        "displacement_distance", displacement_distance,
        "symmetry_tol", symmetry_tol,
        "masses", masses,
        "calc_raman", calc_raman,
        "calc_irreps", calc_irreps,
        "polarization_axes", polarization_axes,
        "calc_bands", calc_bands,
        "calc_displacements", calc_displacements,
        "calc_force_sets", calc_force_sets,
        "write_all_mode_anim", write_all_mode_anim,
        "write_all_mode_gplot", write_all_mode_gplot,
        "supercell_dim", supercell_dim,
        "qpoints", qpoints
    );

    return true;
}

bool sp2::phonopy::phonopy_settings_t::deserialize(const Json::Value &input)
{
    // TODO: make this not a hack, add to serialize()
    auto pol_str = check_pol_axes(input["polarization_axes"]);
    if (pol_str == "backscatter_avg")
    {
        calc_raman_backscatter_avg = true;
    }

    io::deserialize_basic(input,
        "n_samples", n_samples,
        "minimize", min_set,
        "displacement_distance", displacement_distance,
        "symmetry_tol", symmetry_tol,
        "masses", masses,
        "calc_raman", calc_raman,
        "calc_irreps", calc_irreps,
        "polarization_axes", polarization_axes,
        "calc_bands", calc_bands,
        "calc_displacements", calc_displacements,
        "calc_force_sets", calc_force_sets,
        "write_all_mode_anim", write_all_mode_anim,
        "write_all_mode_gplot", write_all_mode_gplot,
        "supercell_dim", supercell_dim,
        "qpoints", qpoints
    );

    return true;
}

bool sp2::phonopy::qpoint_t::deserialize(const Json::Value &input)
{
    if (!input.isArray() || input.size() != 4)
        return false;

    io::get_json_as_type(input[0], label);
    io::get_json_as_type(input[1], x);
    io::get_json_as_type(input[2], y);
    io::get_json_as_type(input[3], z);

    return true;
}

bool sp2::phonopy::qpoint_t::serialize(Json::Value &output) const
{
    output = Json::arrayValue;

    output.append(label);
    output.append(x);
    output.append(y);
    output.append(z);

    return true;
}
