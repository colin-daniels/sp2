#include "program_args_t.hpp"

#include <cxxopts.hpp>

class arg_parser_t
{
    cxxopts::Options options;

    bool display_help = false,
        display_defaults = false,
        run_tests = false;

    bool quiet = false,
        verbose = false;

    std::string config_filename = "";

public:
    arg_parser_t();

    void parse(int &argc, char **&argv);
    std::string get_help_string() { return options.help(); }

    sp2::program_mode get_program_mode();
    sp2::log_level get_log_level();
    std::string get_config_filename() { return config_filename; }
};

arg_parser_t::arg_parser_t() :
    options("sp2", " - Commandline physics program package.")
{
    options.positional_help("[Config File]");

    options.add_options()
        ("h,help",    "Display this information.",
            cxxopts::value<bool>(display_help))
#ifdef SP2_ENABLE_TESTS
        ("t,test",    "Run tests via gtest.",
            cxxopts::value<bool>(run_tests))
#endif // SP2_ENABLE_TESTS
        ("q,quiet",   "Suppress stdout.",
            cxxopts::value<bool>(quiet))
        ("v,verbose", "Display additional information during run.",
            cxxopts::value<bool>(verbose))
        ("c,config",  "Input configuration file.",
            cxxopts::value<std::string>(config_filename))
        ("g,generate-defaults",  "Output default configuration to stdout.",
            cxxopts::value<bool>(display_defaults));

    options.parse_positional("config");
}

void arg_parser_t::parse(int &argc, char **&argv)
{
    options.parse(argc, argv);

    if (argc > 1)
        throw cxxopts::OptionException("Too many positional arguments.");
}

sp2::program_mode arg_parser_t::get_program_mode()
{
    // can only deal with one mode at a time
    int modes_specified = (run_tests + display_help + display_defaults);
    if (modes_specified > 1)
        throw cxxopts::OptionException(
            "Conflicting options specified (e.g. --help and --test).");

    if (run_tests)
        return sp2::program_mode::run_tests;
    if (display_help)
        return sp2::program_mode::display_help;
    if (display_defaults)
        return sp2::program_mode::display_defaults;

    return sp2::program_mode::run_normal;
}

sp2::log_level arg_parser_t::get_log_level()
{
    if (quiet && verbose)
        throw cxxopts::OptionException(
            "Conflicting options --quiet and --verbose.");

    if (quiet)
        return sp2::log_level::error;
    if (verbose)
        return sp2::log_level::info;

    // default minimum output level is warning
    return sp2::log_level::warning;
}

bool find_test_arg(int argc, char **argv)
{
    for (int i = 1; i < argc; ++i)
    {
        std::string arg{argv[i]};
        if (arg == "--test" || arg == "-t")
        {
            return true;
        }
        else if (arg == "--")
            return false;
    }

    return false;
}

sp2::program_args_t sp2::parse_args(int &argc, char **&argv)
{
    arg_parser_t parser;

    program_args_t processed_args;
    processed_args.help_string = parser.get_help_string();

#ifdef SP2_ENABLE_TESTS
    // first try to find -t, or --test since gtest options are not compatible
    // with how we use cxxopts
    if (find_test_arg(argc, argv))
    {
        processed_args.mode = program_mode::run_tests;

        processed_args.argc = argc;
        processed_args.argv = argv;
        return processed_args;
    }
#endif // SP2_ENABLE_TESTS

    try
    {
        // note: modifies argc, argv
        parser.parse(argc, argv);

        processed_args.mode = parser.get_program_mode();
        processed_args.output_level = parser.get_log_level();
        processed_args.config_filename = parser.get_config_filename();

        processed_args.argc = argc;
        processed_args.argv = argv;
    }
    catch (const cxxopts::OptionException& e)
    {
        // rethrow after adding help message to exception,
        // caller will deal with it
        throw std::runtime_error(
            e.what() + ("\n\n" + processed_args.help_string)
        );
    }

    return processed_args;
}
