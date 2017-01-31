//
// Created by cc on 11/25/16.
//

#ifndef SP2_RUN_SETTINGS_T_HPP
#define SP2_RUN_SETTINGS_T_HPP

#include "atac/settings_t.hpp"
#ifdef SP2_ENABLE_PHONOPY
#include "phonopy/phonopy_settings.hpp"
#endif // SP2_ENABLE_PHONOPY
#include "symm/symm_settings_t.hpp"
#include "common/io/structure.hpp"
#include "common/structure_t.hpp"
#include "common/enums.hpp"
#include "common/json/json_serializable_t.hpp"

namespace sp2 {

struct run_settings_t : public io::json_serializable_t
{
    run_type mode = run_type::NONE;
    structure_t structure;

    bool add_hydrogen = false;
    std::string log_filename = "progress.log";

    symm::symm_settings_t symm_settings;
    atac::atac_settings_t atac_settings;
#ifdef SP2_ENABLE_PHONOPY
    phonopy::phonopy_settings_t phonopy_settings;
#endif // SP2_ENABLE_PHONOPY
    minimize::acgsd_settings_t minimize_settings;

    bool serialize(Json::Value &output) const;
    bool deserialize(const Json::Value &input);
};

} // namespace sp2


#endif // SP2_RUN_SETTINGS_T_HPP
