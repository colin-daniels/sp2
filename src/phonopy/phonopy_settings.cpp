//
// Created by cc on 11/26/16.
//

#include "phonopy_settings.hpp"
#include "json/json.h"
#include "common/json/json.hpp"
#include <utility>

using namespace std;

bool sp2::phonopy::phonopy_settings_t::serialize(Json::Value &output) const
{
    io::serialize_basic(output,
        "n_samples", n_samples,
        "minimize", min_set,
        "displacement_distance", displacement_distance,
        "log_filename", log_filename,
        "supercell_dim", supercell_dim,
        "qpoints", qpoints
    );

    return true;
}

bool sp2::phonopy::phonopy_settings_t::deserialize(const Json::Value &input)
{
    io::deserialize_basic(input,
        "n_samples", n_samples,
        "minimize", min_set,
        "displacement_distance", displacement_distance,
        "log_filename", log_filename,
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
