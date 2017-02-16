#ifndef SP2_RELAXATION_SETTINGS_HPP
#define SP2_RELAXATION_SETTINGS_HPP

#include <common/minimize/settings.hpp>
#include <common/json/json_serializable_t.hpp>

namespace sp2 {

struct relaxation_settings_t : public io::json_serializable_t
{
    std::string output_file;
    minimize::acgsd_settings_t minimize_settings;

    bool serialize(Json::Value &output) const;
    bool deserialize(const Json::Value &input);
};

} // namespace sp2

#endif // SP2_RELAXATION_SETTINGS_HPP
