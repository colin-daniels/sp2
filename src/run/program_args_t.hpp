#ifndef SP2_PROGRAM_ARGS_T_HPP
#define SP2_PROGRAM_ARGS_T_HPP

#include <string>
//#include <common/env/env_t.hpp>

namespace sp2 {

enum class log_level
{
    info,
    warning,
    error
};

enum class program_mode
{
    none,
    display_help,
    display_defaults,
    run_tests,
    run_normal
};

struct program_args_t
{
    int argc = 0;
    char **argv = nullptr;

    program_mode mode = program_mode::none;
    log_level output_level = log_level::warning;

    std::string config_filename = "";
    std::string help_string = "";
};

/// parse program arguments from input
program_args_t parse_args(int &argc, char **&argv);

} // namespace sp2

#endif // SP2_PROGRAM_ARGS_T_HPP
