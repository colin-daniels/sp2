#include "src/run/run_types.hpp"
#include "src/common/minimize/minimize.hpp"

#include "src/lammps/lammps_interface.hpp"
#include "src/airebo/system_control_t.hpp"
#include "src/phos/ephos.hpp"

int sp2::run_relaxation(const run_settings_t &config, MPI_Comm)
{
    if (config.relaxation_settings.output_file.empty())
    {
        std::cerr << "No output file name given for structural relaxation."
                  << std::endl;

        return EXIT_FAILURE;
    }

    auto run_relax = [&](auto&& sys) {
        // minimize
        sys.set_position(
            minimize::acgsd(sys.get_diff_fn(), sys.get_position(),
                config.relaxation_settings.minimize_settings)
        );

        // write output
        if (!io::write_structure(config.relaxation_settings.output_file,
                sys.get_structure()))
        {
            return EXIT_FAILURE;
        }

        return EXIT_SUCCESS;
    };

    switch (config.potential)
    {
    case potential_type::LAMMPS:
        return run_relax(lammps::system_control_t(
            config.structure, config.lammps_settings));
    case potential_type::REBO:
        return run_relax(airebo::system_control_t(config.structure));
    case potential_type::PHOSPHORENE:
        return run_relax(phos::phosphorene_sys_t(config.structure));
    default:
        return EXIT_FAILURE;
    }
}
