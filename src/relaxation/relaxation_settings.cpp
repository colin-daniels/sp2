#include <common/json/json.hpp>
#include "relaxation_settings.hpp"

void sp2::relaxation_settings_t::serialize(Json::Value &output) const
{
    io::serialize_basic(output,
        "output_file", output_file,
        "minimize", minimize_settings
    );
}

bool sp2::relaxation_settings_t::deserialize(const Json::Value &input)
{
    return io::deserialize_basic(input,
        "output_file", output_file,
        "minimize", minimize_settings
    );
}
