
#warning FIXME need to similarly guard out use-sites of bindings.h so that this can compile without the flag
#ifdef SP2_ENABLE_PYTHON
#include "Python.h"
#include "common/python/bindings.h"
#endif // SP2_ENABLE_PYTHON

#include <iostream>
#include <functional>
#include <map>

#include "run/run_types.hpp"
#include "common/json/json.hpp"

#include <boost/mpi.hpp>
#include <common/util/templates.hpp>

#ifdef SP2_ENABLE_TESTS
#include <gtest/gtest.h>
#endif // SP2_ENABLE_TESTS

using namespace std;
using namespace sp2;

/// normal operation, switches run type based on input config
int run_normal(string config_filename, MPI_Comm comm);

int main(int argc, char *argv[])
{
    // setup program environment
    boost::mpi::environment env;

    #ifdef SP2_ENABLE_PYTHON
    sp2::python::initialize(argv[0]);

    // FIXME don't use atexit
    auto guard = scope_guard([&] { sp2::python::finalize(); });
    #endif // SP2_ENABLE_PYTHON

////////////////////////////////////////////////////////////////////////////////
// check command line flags for special options that bypass normal operation  //
////////////////////////////////////////////////////////////////////////////////

    // TODO: parse args better, verbose logging
    map<string, function<int(void)>> run_special;

    run_special["--generate-defaults"] = []{
        // note: declared in run/run_types.hpp
        return generate_defaults("config.json");
    };

#ifdef SP2_ENABLE_TESTS
    run_special["--test"] = [&]{
        ::testing::InitGoogleTest(&argc, argv);
        return RUN_ALL_TESTS();
    };
#endif // SP2_ENABLE_TESTS

    for (int i = 1; i < argc; ++i)
        if (run_special.count(argv[i]))
            return run_special[argv[i]]();

////////////////////////////////////////////////////////////////////////////////
// normal operation                                                           //
////////////////////////////////////////////////////////////////////////////////

    // get config filename (if specified)
    string config_filename = "config.json";
    if (argc > 1)
        config_filename = string(argv[1]);

    // start a normal run
    return run_normal(config_filename, MPI_COMM_WORLD);
}

int run_normal(string config_filename, MPI_Comm comm)
{
    // read the json configuration file
    Json::Value config;
    if (!io::read_json_file(config, config_filename))
    {
        cerr << "Failed to load json configuration file, run with the flag "
             << "--generate-defaults to generate a default configuration file "
             << "with the name 'config.json' in the current directory."
             << endl;
        return EXIT_FAILURE;
    }

////////////////////////////////////////////////////////////////////////////////
// determine the type of run and execute it                                   //
////////////////////////////////////////////////////////////////////////////////

    // try to read all the settings from the config
    run_settings_t settings;
    if (!settings.deserialize(config))
    {
        cerr << "Failed to deserialize json configuration file \""
             << config_filename << "\"" << endl;
        return EXIT_FAILURE;
    }


    // execute the run (functions are declared in run/run_types.hpp)
    switch (settings.mode)
    {
    case run_type::RELAX:
        // just minimize/relax the structure and nothing else
        return run_relaxation(settings, comm);
    case run_type::ATAC:
        // run ATAC
//        return run_atac(settings, comm);
    case run_type::SYMM:
        // symmetry structure search
        return run_symm(settings, comm);
#ifdef SP2_ENABLE_PHONOPY
    case run_type::PHONOPY:
        // phonopy
        return run_phonopy(settings, comm);
#endif // SP2_ENABLE_PHONOPY
    default:
        return EXIT_FAILURE;
    }
}
