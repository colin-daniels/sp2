#ifndef SP2_ATAC_SETTINGS_T_HPP
#define SP2_ATAC_SETTINGS_T_HPP

#include "common/json/json_serializable_t.hpp"
#include "common/minimize/settings.hpp"

namespace sp2 {
namespace atac {

struct atac_settings_t : public io::json_serializable_t
{
    /// number of iterations to run ATAC, 0 for no limit
    int iteration_limit;

    double temperature, ///< simulation temperature in Kelvin
        rate_stw,       ///< Stone-Thrower-Wales probability [0:1]
        rate_divacancy, ///< two atom vacancy (divacancy) probability [0:1]
        rate_add_dimer, ///< two adatom (add dimer) probability [0:1]
        rate_cyclo;     ///< cycloaddition (sp2->sp3->sp2) probability [0:1]

    /// minimum neighbor distance for initial bond graph generation
    double bond_min;
    /// maximum neighbor distance for initial bond graph generation
    double bond_max;

    double penalty_scale,   ///< penalty potential multiplier (eV/A^2)
        penalty_cutoff;     ///< penalty potential distance cutoff (A)

    // std::string output_filename;

    int stretch_interval;
    double stretch_amount[3];

    /// ATAC specific minimization settings
    minimize::acgsd_settings_t min_set;

    atac_settings_t();

    void serialize(Json::Value &output) const;

    bool deserialize(const Json::Value &input);
};

} // namespace atac
} // namespace sp2
#endif // SP2_ATAC_SETTINGS_T_HPP
