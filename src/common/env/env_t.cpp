#include "env_t.hpp"

#ifdef SP2_ENABLE_TESTS
#include <gtest/gtest.h>

#endif // SP2_ENABLE_TESTS
#include <run/run_settings_t.hpp>

sp2::env_t::env_t(const boost::mpi::communicator &comm_in,
    const logger_t &logger_in, randutils::mt19937_rng  &&rng_in) :
        comm(comm_in), logger(logger_in), rng(rng_in) {}


// TODO: enable/disable MPI
namespace mpi = boost::mpi;

int run(int argc, char *argv[]);
int run(int argc, char *argv[], MPI_Comm comm_in);
int run(sp2::program_args_t pargs, mpi::communicator comm);

int display_help(sp2::program_args_t pargs);
int display_defaults();

int run_tests(sp2::program_args_t pargs);
int run_normal(sp2::program_args_t pargs, mpi::communicator comm);

void verify_mpi_environment();

int run(int argc, char *argv[])
{
    return run(argc, argv, MPI_COMM_WORLD);
}

int run(int argc, char *argv[], MPI_Comm comm_in)
{
    // Initializes the mpi environment (i.e. via MPI_Init or similar), when
    // destroyed, calls MPI_Finalize/Abort if a result of an exception.
    //
    // Note: does nothing if the mpi environment has already been initialized
    auto mpi_env = std::make_unique<mpi::environment>(argc, argv,
        // Only main thread will do MPI calls (equiv to MPI_THREAD_FUNNELED)
        mpi::threading::level::funneled);

    // In the case that the environment has already been set up, we need to
    // verify that it's what we expect
    verify_mpi_environment();

    // attach to the input communicator (relies on calling code to
    // manage/free the input communicator)
    mpi::communicator comm(comm_in, mpi::comm_create_kind::comm_attach);
    if (!comm)
        throw std::runtime_error("Invalid input MPI_Comm.");


    // Read program arguments
    sp2::program_args_t pargs;
    try {
        pargs = sp2::parse_args(argc, argv);
    } catch (const std::exception &e) {
        std::cerr << "Error parsing options: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    return run(pargs, comm);
}

void verify_mpi_environment()
{
    if (mpi::environment::finalized())
        throw std::runtime_error("MPI has already been finalized.");

    if (!mpi::environment::initialized())
        throw std::runtime_error("MPI has not been initialized.");

    // our code is written such that it assumes the threading level is
    // funneled, so we need to check it
    switch (mpi::environment::thread_level())
    {
    case mpi::threading::level::single:
        throw std::runtime_error("MPI_THREAD_SINGLE is not supported.");

        // funneled is equivalent to serialized in our case, so its O.K.
    case mpi::threading::level::funneled:
    case mpi::threading::level::serialized:
        break;

    case mpi::threading::level::multiple:
        // known to be buggy, at least for OpenMPI
        throw std::runtime_error("MPI_THREAD_MULTIPLE is not supported.");
    }
}

int run(sp2::program_args_t pargs, mpi::communicator comm)
{
    // Switch depending on the "program mode",
    switch (pargs.mode)
    {
    case sp2::program_mode::display_help:     return display_help(pargs);
    case sp2::program_mode::display_defaults: return display_defaults();
    case sp2::program_mode::run_tests:        return run_tests(pargs);
    case sp2::program_mode::run_normal:       return run_normal(pargs, comm);

    case sp2::program_mode::none:
        break;
    }

    std::cerr << "Invalid or unspecified program mode.\n";
    return EXIT_FAILURE;
}


int display_help(sp2::program_args_t pargs)
{
    std::cout << pargs.help_string << std::endl;
    return EXIT_SUCCESS;
}

int display_defaults()
{
    Json::Value output;

    // the default constructed values of the fields in the settings structures
    // are the actual defaults for the program
    sp2::run_settings_t().serialize(output);

    std::string content;
    if (!sp2::io::write_json(output, content))
        return EXIT_FAILURE;

    std::cout << content << std::endl;
    return EXIT_SUCCESS;
}

int run_tests(sp2::program_args_t pargs)
{
#ifdef SP2_ENABLE_TESTS
    ::testing::InitGoogleTest(&pargs.argc, pargs.argv);
    return RUN_ALL_TESTS();
#else
    (void)(pargs);
    std::cerr << "Cannot run tests, recompile with them enabled.\n";
    return EXIT_FAILURE;
#endif // SP2_ENABLE_TESTS
}

sp2::run_settings_t read_config(sp2::program_args_t pargs)
{
    // read the json configuration file
    Json::Value config;
    if (!sp2::io::read_json_file(config, pargs.config_filename))
    {
        throw std::runtime_error(
            "Failed to load json configuration file, run with the flag"
            " --generate-defaults to output the default configuration."
        );
    }

    // try to read all the settings from the config
    sp2::run_settings_t settings;
    if (!settings.deserialize(config))
        throw std::runtime_error("Failed to deserialize json configuration file"
            "\"" + pargs.config_filename + "\"");

    return settings;
}

int run_normal(sp2::program_args_t pargs, mpi::communicator)
{
    // read configuration file
    sp2::run_settings_t settings;
    try {
        settings = read_config(pargs);
    } catch (const std::runtime_error &e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    // TODO: yaml input, command-like structure?

    // TODO: setup environment, const settings in env?

    // TODO: switch on run type

    return EXIT_SUCCESS;
}


#ifdef SP2_ENABLE_TESTS

TEST(env, all) {
//    boost::mpi::environment env(argc, argv,
//        boost::mpi::threading::level::funneled);

//    auto env = sp2::construct_thread_environments(3);
//
//    randutils::seed_seq_fe256 sseq{
//        0x97
//    };
//
//    auto env2 = sp2::construct_thread_environments(3, {}, {}, std::move(sseq));
}

#endif // SP2_ENABLE_TESTS
