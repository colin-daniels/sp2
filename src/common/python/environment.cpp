#include <Python.h>

#include "environment.hpp"
#include "fake_modules.hpp"
#include "include_numpy.hpp"
#include "ext_modules.hpp"
#include "error.hpp"
#include "py_scoped_t.hpp"

#include <iostream>

//-----------------------------

using namespace std;
using namespace sp2::python;
namespace sp2 {
namespace python {

// anonymous namespace for private things
namespace {

void extend_sys_path_single(const char *dir)
{

    // python literal equivalent:
    //
    // if d in set(sys.path):
    //     sys.path.remove(d)
    // sys.path.insert(0, d)

    auto py_dir = scope(PyUnicode_DecodeFSDefault(dir));
    throw_on_py_err("add_to_sys_path: error decoding filesystem path");

    auto sys_path = scope_dup(PySys_GetObject((char *) "path"));
    throw_on_py_err("add_to_sys_path: error reading sys.path");

    // make a set because I can't find where the C API exposes list.__contains__... >_>
    auto sys_path_set = scope(PySet_New(sys_path.raw()));
    throw_on_py_err("add_to_sys_path: error making set from sys.path");

    // remove an existing entry
    int already_there = PySet_Contains(sys_path_set.raw(), py_dir.raw());
    throw_on_py_err("add_to_sys_path: error checking set membership");
    switch (already_there)
    {
    case 1:
        PyObject_CallMethod(sys_path.raw(), "remove", "o", py_dir.raw());
        throw_on_py_err("add_to_sys_path: error removing existing item");
        break;
    case 0:
        break;
    default:
        throw logic_error("unexpected result from __contains__");
    }

    // prepend
    PyList_Insert(sys_path.raw(), 0, py_dir.raw());
    throw_on_py_err("add_to_sys_path: error inserting item");
}

} // anonymous namespace

/* --------------------------------------------------------------------- */
// Public API

void environment::ensure_run_once()
{
    static bool has_run = false;
    if (has_run)
        throw logic_error("Initialized py interpreter more than once!");
    has_run = true;
}

environment::environment(const char *prog)
{
    ensure_run_once();

    if (prog)
    {
        py_allocated_program = Py_DecodeLocale(prog, NULL);
        if (py_allocated_program)
        {
            Py_SetProgramName(py_allocated_program);
        }
        else
        {
            throw runtime_error(
                "Warning: Could not decode program name for python bindings");
        }
    }

    ext_modules::pre_py_initialize();
    Py_Initialize();
    ext_modules::post_py_initialize();

    initialize_numpy();
    throw_on_py_err("error initializing numpy");

    fake_modules::initialize();
}

int environment::finalize()
{
    fake_modules::finalize();
    if (int code = Py_FinalizeEx())
        return code;
    if (py_allocated_program)
        PyMem_RawFree(py_allocated_program);
    return 0;
}

environment::~environment()
{
    try
    {
        if (finalize())
        {
            std::cerr << "WARNING: Unknown errors occurred in Python while"
                      << " cleaning up the python interpreter." << std::endl;
        }
    }
    catch (const std::exception &e)
    {
        // don't throw from destructor
        std::cerr << "WARNING: An uncaught exception was thrown while cleaning"
                  << " up the python interpreter: " << e.what() << std::endl;
    }
    catch (...)
    {
        // I said, /don't throw from destructor/!!!
        std::cerr << "WARNING: Something ugly and unrecognizable was thrown"
                  << " our way while cleaning up the python interpreter"
                  << std::endl;
    }
}

void extend_sys_path(vector<string> dirs)
{
    // reverse order so that the prepended results appear in the requested order
    for (auto it = dirs.rbegin(); it != dirs.rend(); it++)
        extend_sys_path_single(it->c_str());
}

} // namespace sp2
} // namespace python
