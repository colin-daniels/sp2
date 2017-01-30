#include "run/run_types.hpp"

#include "common/io/structure.hpp"
#include "atac/settings_t.hpp"
#include "common/minimize/settings.hpp"
#include "lammps/settings_t.hpp"
#include "common/json/json.hpp"

using namespace std;
using namespace sp2;

int sp2::generate_defaults(std::string filename)
{
    Json::Value output;

    // the default constructed values of the fields in the settings structures
    // are the actual defaults for the program
    run_settings_t().serialize(output);

    if (!io::write_json_file(output, filename))
        return EXIT_FAILURE;

    return EXIT_SUCCESS;
}
