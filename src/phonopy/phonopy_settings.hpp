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

struct phonopy_metro_settings_t : public io::json_serializable_t
{
    bool enabled = false;

    // Directories to be prepended to sys.path, where python modules may be found.
    //
    // The default prepends an entry of "" (the current directory), just like
    // the python interpreter itself typically does behind the scenes.
    std::vector<std::string> python_sys_path = {""};
    std::string python_module = "mutate";
    std::string python_function = "mutate";
    minimize::metropolis_settings_t settings;

    virtual bool serialize(Json::Value &output) const;
    virtual bool deserialize(const Json::Value &input);
};

struct phonopy_settings_t : public io::json_serializable_t
{
    int n_samples = 250;
    int supercell_dim[3] = {1, 1, 1};
    double displacement_distance = 0.01,
        symmetry_tol = 1e-5;

    bool calc_raman = false,
        calc_raman_backscatter_avg = false,
        calc_irreps = false,
        write_all_mode_anim = false,
        write_all_mode_gplot = false;

    std::vector<double> masses;

    // TODO: serialize
    double raman_active_cutoff = 1e-3;
    int polarization_axes[2] = {0, 0};

    bool calc_displacements = true;
    bool calc_force_sets = true;
    bool calc_bands = true;

    std::vector<qpoint_t> qpoints = {
        qpoint_t("\u0393", 0, 0, 0),
        qpoint_t("X", 0.5, 0, 0)
    };

    minimize::acgsd_settings_t min_set;
    phonopy_metro_settings_t metro_set;

    virtual bool serialize(Json::Value &output) const;
    virtual bool deserialize(const Json::Value &input);
};

} // namespace phonopy
} // namespace sp2

#endif // SP2_PHONOPY_SETTINGS_T_HPP
