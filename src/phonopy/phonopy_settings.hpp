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

struct phonopy_metro_funcs_t : public io::json_serializable_t {
    bool advanced = true;

    /// Name for an all-in-one mutation function.
    /// Only meaningful if 'advanced == false'.
    ///
    /// This should return a value which can be recognized
    ///  as a 'structural_mutation_t'.
    std::string mutate = "mutate";

    /// Name for a mutation generation function.
    /// Only meaningful if 'advanced == true'.
    ///
    /// If defined, then it should accept the same arguments
    /// as 'mutate', but can return any arbitrary python object.
    /// The other callbacks will decide how to interpret this object.
    std::string generate = "generate";

    /// Name for a mutation-applying function.
    /// Only meaningful if 'advanced == true'.
    ///
    /// If defined, then it should accept:
    ///  (1) the output of generate, as a positional argument, followed by
    ///  (2) the standard kw arguments provided to 'mutate'
    /// It should produce a 'structural_mutation_t'
    ///
    /// If omitted, will search for a function named 'apply',
    /// and then if that is not found, the following definition is assumed:
    ///
    /// def apply(mutation, **kw):
    ///     return mutation  # assume mutation is 'structural_mutation_t'
    std::string apply;

    /// Name for a 'mutation accepted' callback.
    /// Only meaningful if 'advanced == true'.
    ///
    /// If defined, then it should accept:
    ///  (1) the output of generate, as a positional argument, followed by
    ///  (2) the standard kw arguments provided to 'mutate'
    /// It should produce 'None'.
    ///
    /// This function will be called when a mutation is accepted by the
    /// metropolis algorithm.
    ///
    /// If omitted, will search for a function named 'apply',
    /// and then if that is not found, the following definition is assumed:
    ///
    /// def accept(mutation, **kw):
    ///     pass
    std::string accept;

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
    phonopy_metro_funcs_t python_functions;
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
