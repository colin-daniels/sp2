//
// Created by cc on 11/26/16.
//

#include "phonopy_settings.hpp"
#include "common/json/json.hpp"
#include <utility>

using namespace std;

Json::Value serialize_polarization(const sp2::phonopy::phonopy_settings_t &obj)
{
    if (obj.calc_raman_backscatter_avg) {
        return "backscatter_avg";
    } else {
        Json::Value output = Json::arrayValue;
        output.append(obj.polarization_axes[0]);
        output.append(obj.polarization_axes[1]);
        return output;
    }
}

bool deserialize_polarization(
    sp2::phonopy::phonopy_settings_t &obj,
    const Json::Value &input)
{
    switch (input.type()) {
        case Json::stringValue:
            if (input.asString() == "backscatter_avg") {
                obj.calc_raman_backscatter_avg = true;
                return true;
            }
            return false;

        case Json::arrayValue:
            if (input.size() != 2)
                return false;

            return sp2::io::get_json_as_type(input[0], obj.polarization_axes[0])
                && sp2::io::get_json_as_type(input[1], obj.polarization_axes[1]);

        default:
            return false;
    }
}

void sp2::phonopy::phonopy_settings_t::serialize(Json::Value &output) const
{
    Json::Value pol = serialize_polarization(*this);

    io::serialize_basic(output,
        "n_samples", n_samples,
        "minimize", min_set,
        "displacement_distance", displacement_distance,
        "symmetry_tol", symmetry_tol,
        "masses", masses,
        "calc_raman", calc_raman,
        "calc_irreps", calc_irreps,
        "compute_force_constants", compute_force_constants,
        "polarization_axes", pol,
        "do_minimization", do_minimization,
        "calc_bands", calc_bands,
        "calc_displacements", calc_displacements,
        "calc_force_sets", calc_force_sets,
        "metropolis", metro_set,
        "write_all_mode_anim", write_all_mode_anim,
        "write_all_mode_gplot", write_all_mode_gplot,
        "supercell_dim", supercell_dim,
        "write_force", write_force,
        "qpoints", qpoints
    );
}

bool sp2::phonopy::phonopy_settings_t::deserialize(const Json::Value &input)
{
    Json::Value pol;

    if (!deserialize_polarization(*this, pol))
        return false;

    return io::deserialize_basic(input,
        "n_samples", n_samples,
        "minimize", min_set,
        "displacement_distance", displacement_distance,
        "symmetry_tol", symmetry_tol,
        "masses", masses,
        "calc_raman", calc_raman,
        "calc_irreps", calc_irreps,
        "polarization_axes", pol,
        "compute_force_constants", compute_force_constants,
        "do_minimization", do_minimization,
        "calc_bands", calc_bands,
        "calc_displacements", calc_displacements,
        "calc_force_sets", calc_force_sets,
        "metropolis", metro_set,
        "write_all_mode_anim", write_all_mode_anim,
        "write_all_mode_gplot", write_all_mode_gplot,
        "supercell_dim", supercell_dim,
        "write_force", write_force,
        "qpoints", qpoints
    );
}

bool sp2::phonopy::qpoint_t::deserialize(const Json::Value &input)
{
    if (!input.isArray() || input.size() != 4)
        return false;

    return io::get_json_as_type(input[0], label)
           && io::get_json_as_type(input[1], x)
           && io::get_json_as_type(input[2], y)
           && io::get_json_as_type(input[3], z);
}

void sp2::phonopy::qpoint_t::serialize(Json::Value &output) const
{
    output = Json::arrayValue;

    output.append(label);
    output.append(x);
    output.append(y);
    output.append(z);
}
