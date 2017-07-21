//
// Created by lampam on 7/11/17.
//

// Python.h must be included prior to any other headers.
#include "Python.h"
#include <common/python/bindings.hpp>

#include <common/python/util.hpp>
#include <common/python/conversion.hpp>
#include <common/python/include_numpy.hpp>
#include <common/python/modules/fake_modules.hpp>

//-----------------------------

using namespace std;
using namespace sp2::python;

// namespace for private things
namespace sp2 {
namespace python {
namespace impl {

/* --------------------------------------------------------------------- */
// helper type for python function arguments

// a pair of *args and **kw
typedef pair<py_scoped_t, py_scoped_t> args_t;

args_t args_and_kw(py_scoped_t args, py_scoped_t kw)
{
    return make_pair(move(args), move(kw));
}

args_t just_args(py_scoped_t args)
{
    auto kw = scope(PyDict_New());
    return make_pair(move(args), move(kw));
}

args_t just_kw(py_scoped_t kw)
{
    auto args = scope(PyTuple_New(0));
    return make_pair(move(args), move(kw));
}

/* --------------------------------------------------------------------- */
// helper functions

// Access an attribute of a python object.
//
// Equivalent to 'getattr(obj, name)'.
// Throws an exception if the attribute does not exist.
py_scoped_t getattr(py_scoped_t &o, const char *attr)
{
    auto tmp = scope(PyObject_GetAttrString(o.raw(), attr));
    throw_on_py_err();
    return move(tmp);
}

// Invoke an object's __call__.
py_scoped_t call_callable(py_scoped_t &function, args_t argpair)
{
    auto &args = argpair.first;
    auto &kw = argpair.second;
    if (!(args && kw))
        throw logic_error(
            "tried to call function with invalid or previously used args_t");

    if (!PyArg_ValidateKeywordArguments(kw.raw()))
        throw logic_error("invalid keyword arguments");
    throw_on_py_err(); // errors not covered by the above

    if (!PyCallable_Check(function.raw()))
        throw logic_error("not callable");

    auto retval = scope(PyObject_Call(function.raw(), args.raw(), kw.raw()));
    throw_on_py_err("error calling python callable");
    return retval;
}

// Call a named function in a named module with *args and **kw.
// The module will be automatically imported if it isn't already loaded;
//   in particular, calling this multiple times will not result in repeated
//   calls to the module's __init__.
//
// Equivalent to:
//
//     import mod_name
//     return mod_name.func_name(*args, **kw)
//
py_scoped_t call_module_function(const char *mod_name, const char *func_name,
    args_t args)
{
    py_scoped_t module;
    {
        auto name = scope(PyUnicode_DecodeFSDefault(mod_name));
        throw_on_py_err();

        module = scope(PyImport_Import(name.raw()));
        throw_on_py_err();
    }

    auto func = getattr(module, func_name);
    if (!(PyCallable_Check(func.raw())))
        throw runtime_error(string() + mod_name + "."
                            + func_name + " is not callable");

    return call_callable(func, move(args));
}

sp2::structural_mutation_t call_run_phonopy_mutation_function(
    const char *mod_name, const char *func_name,
    vector<double> carts, double lattice[3][3],
    vector<size_t> sc_to_prim)
{
    // interpret as 3N cartesian coords
    size_t width = 3;
    if (carts.size() % width != 0)
        throw logic_error("vector not divisible by width");

    size_t height = carts.size() / width;

    auto py_carts = to_python_strict(as_ndarray(carts, {height, width}));

    auto c_lattice = vector<double>(&lattice[0][0], &lattice[0][0] + 9);
    auto py_lattice = to_python_strict(as_ndarray(c_lattice, {3, 3}));

    py_scoped_t py_sc_map = [&] {
        py_scoped_t list = to_python_strict(sc_to_prim);
        auto &module = fake_modules::mutation_helper.module;
        auto klass = getattr(module, "supercell_index_mapper");

        auto args = scope(Py_BuildValue("(O)", list.raw()));
        throw_on_py_err("Exception constructing python args tuple.");

        return call_callable(klass, just_args(args.dup()));
    }();

    py_scoped_t kw = [&] {
        auto kw = scope(Py_BuildValue("{sOsOsO}",
            "carts", py_carts.raw(),
            "lattice", py_lattice.raw(),
            "supercell", py_sc_map.raw()
        ));
        throw_on_py_err("Exception constructing python kw dict.");
        return move(kw);
    }();

    auto py_retval = call_module_function(mod_name, func_name,
        just_kw(kw.dup()));

    return from_python_strict<structural_mutation_t>(py_retval,
        "error converting python return value");
}

void extend_sys_path(const char *dir)
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


void extend_sys_path(vector<string> dirs)
{
    // reverse order so that the prepended results appear in the requested order
    for (auto it = dirs.rbegin(); it != dirs.rend(); it++)
        extend_sys_path(it->c_str());
}

void initialize(const char *prog)
{
    if (prog)
    {
        // NOTE: Technically, we are leaking this memory by not setting up any sort of
        //       provision to call 'PyMem_RawFree(program)' after 'Py_FinalizeEx()'.
        //       But this is no big deal since this function should only be called once.
        wchar_t *program = Py_DecodeLocale(prog, NULL);
        if (program)
        {
            Py_SetProgramName(program);
        }
        else
        {
            throw runtime_error(
                "Warning: Could not decode program name for python bindings");
        }
    }

    Py_Initialize();

    initialize_numpy();
    throw_on_py_err("error initializing numpy");

    initialize_fake_modules();
}

int finalize()
{
    finalize_fake_modules();
    return Py_FinalizeEx();
}

} // namespace impl
} // namespace python
} // namespace sp2

/* --------------------------------------------------------------------- */
// Public API

void sp2::python::initialize(const char *prog)
{
    return impl::initialize(prog);
}

int sp2::python::finalize()
{
    return impl::finalize();
}

void sp2::python::extend_sys_path(vector<string> dirs)
{
    return impl::extend_sys_path(dirs);
}

sp2::structural_mutation_t sp2::python::call_run_phonopy_mutation_function(
    const char *mod_name, const char *func_name,
    vector<double> carts, double lattice[3][3],
    vector<size_t> sc_to_prim)
{
    return impl::call_run_phonopy_mutation_function(mod_name, func_name, carts,
        lattice, sc_to_prim);
}

