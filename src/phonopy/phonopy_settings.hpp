//
// Created by cc on 11/26/16.
//

#ifndef SP2_PHONOPY_SETTINGS_T_HPP
#define SP2_PHONOPY_SETTINGS_T_HPP

#include <common/json/json_serializable_t.hpp>
#include <vector>
#include <string>
#include <common/minimize/settings.hpp>

namespace sp2 {
namespace phonopy {

struct qpoint_t : public io::json_serializable_t
{
    std::string label = "G";
    double x = 0,
        y = 0,
        z = 0;

    qpoint_t() = default;
    qpoint_t(std::string label_in, double x_in, double y_in, double z_in) :
        label(label_in), x(x_in), y(y_in), z(z_in) {}

    virtual bool serialize(Json::Value &output) const;
    virtual bool deserialize(const Json::Value &input);
};

struct phonopy_settings_t : public io::json_serializable_t
{
    int n_samples = 250;
    int supercell_dim[3] = {1, 1, 1};
    double displacement_distance = 0.01;

    std::vector<qpoint_t> qpoints = {
        qpoint_t("\u0393", 0, 0, 0),
        qpoint_t("X", 0.5, 0, 0)
    };

    std::string log_filename = "progress.log";

    minimize::acgsd_settings_t min_set;

    virtual bool serialize(Json::Value &output) const;
    virtual bool deserialize(const Json::Value &input);
};

} // namespace phonopy
} // namespace sp2

#endif // SP2_PHONOPY_SETTINGS_T_HPP
