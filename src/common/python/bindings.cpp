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
namespace sp2 {
namespace python {

py_scoped_t& sp2::python::py_opaque_t::inner()
{
    if (!_impl) {
        throw std::logic_error("detected py_opaque_t with a NULL unique_ptr!"
            " The inner py_scoped_t should be NULL instead.");
    }

    return *_impl;
}

const py_scoped_t& sp2::python::py_opaque_t::inner() const
{
    if (!_impl) {
        throw std::logic_error("detected py_opaque_t with a NULL unique_ptr!"
            " The inner py_scoped_t should be NULL instead.");
    }

    return *_impl;
}

//
py_opaque_t::py_opaque_t() = default;
py_opaque_t::~py_opaque_t() = default;
py_opaque_t::py_opaque_t(py_opaque_t &&) = default;
py_opaque_t& py_opaque_t::operator=(py_opaque_t &&) = default;
py_opaque_t::py_opaque_t(unique_ptr<py_opaque_t::impl_t> &&impl)
    : _impl(move(impl))
{
    if (!_impl) {
        throw std::logic_error("cannot construct py_opaque_t from "
            "NULL unique_ptr. (the py_scoped_t should be null instead)");
    }
}

// anonymous namespace for private things
namespace {

// less typing to construct 'py_opaque_t's

py_opaque_t opaque(py_scoped_t &&scoped) {
    typedef py_opaque_t::impl_t impl_t;
    return py_opaque_t{make_unique<impl_t>(move(scoped))};
}

py_opaque_t opaque(py_scoped_t &scoped) {
    return opaque(scoped.dup());
}



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

/// Invoke an object's __call__.
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

    cerr << "DEBUG!!" << endl;
    cerr << repr(function) << endl;
    cerr << "ARG: " << repr(args) << endl;
    cerr << "KWS: " << repr(kw) << endl;
    auto retval = scope(PyObject_Call(function.raw(), args.raw(), kw.raw()));
    throw_on_py_err("error calling python callable");
    return retval;
}

// Overload for convenience when there are both kw and args.
// (having only one would be ambiguous)
py_scoped_t call_callable(py_scoped_t &function, py_scoped_t &args,
    py_scoped_t &kw)
{
    auto argpair = args_and_kw(args.dup(), kw.dup());
    return call_callable(function, move(argpair));
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
    auto func = lookup_required(mod_name, func_name).inner().move();

    if (!(PyCallable_Check(func.raw())))
        throw runtime_error(string() + mod_name + "."
                            + func_name + " is not callable");

    return call_callable(func, move(args));
}

py_scoped_t make_run_phonopy_param_pack(
    vector<double> carts, const double lattice[3][3],
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

        cerr << "WOAH LIST: " << repr(list) << endl;
        cerr << "WOAH ARGS: " << repr(args) << endl;

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

    return kw;
}

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

py_scoped_t import_named_module(const char *mod_name) {
    auto name = scope(PyUnicode_DecodeFSDefault(mod_name));
    throw_on_py_err();

    auto module = scope(PyImport_Import(name.raw()));
    throw_on_py_err();

    return module;
}

} // anonymous namespace

/* --------------------------------------------------------------------- */
// Public API

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

void extend_sys_path(vector<string> dirs)
{
    // reverse order so that the prepended results appear in the requested order
    for (auto it = dirs.rbegin(); it != dirs.rend(); it++)
        extend_sys_path_single(it->c_str());
}

py_opaque_t lookup_required(const char *mod_name, const char *attr)
{
    auto module = import_named_module(mod_name);
    return opaque(getattr(module, attr));
}

py_opaque_t lookup_optional(const char *mod_name, const char *attr)
{
    auto module = import_named_module(mod_name);
    return opaque(getattr(module, attr, {}));
}

py_opaque_t run_phonopy::make_param_pack(
    vector<double> carts, const double lattice[3][3],
    vector<size_t> sc_to_prim)
{
    auto py = make_run_phonopy_param_pack(carts, lattice, sc_to_prim);
    return opaque(py);
}

sp2::structural_mutation_t run_phonopy::call_mutate(
    const py_opaque_t &function,
    const py_opaque_t &param_pack)
{
    auto callable = function.inner().dup();
    auto kw = param_pack.inner().dup();
    if (!callable)
        throw logic_error("mutate function must not be NULL");

    auto out = call_callable(callable, just_kw(kw.move()));

    return from_python_strict<structural_mutation_t>(out,
        "error converting python return value");
}

sp2::structural_mutation_t run_phonopy::call_apply(
    const py_opaque_t &function,
    const py_opaque_t &mutation,
    const py_opaque_t &param_pack)
{
    auto callable = function.inner().dup();
    auto py_mutation = mutation.inner().dup();
    auto kw = param_pack.inner().dup();

    py_scoped_t py_cxx_mutation;
    if (callable) {
        auto args = scope(Py_BuildValue("(O)", py_mutation.raw()));
        py_cxx_mutation = call_callable(callable, args, kw);
    } else {
        // default definition: return mutation
        py_cxx_mutation = py_mutation.dup();
    }

    return from_python_strict<structural_mutation_t>(py_cxx_mutation,
        "error converting python return value");
}

void run_phonopy::call_accept(
    const py_opaque_t &function,
    const py_opaque_t &mutation,
    const py_opaque_t &param_pack)
{
    auto callable = function.inner().dup();
    auto py_mutation = mutation.inner().dup();
    auto kw = param_pack.inner().dup();
    if (!callable)
        return; // default definition

    auto args = scope(Py_BuildValue("(O)", py_mutation.raw()));
    call_callable(callable, args, kw);
}

py_opaque_t run_phonopy::call_generate(
    const py_opaque_t &function,
    const py_opaque_t &param_pack)
{
    auto callable = function.inner().dup();
    auto kw = param_pack.inner().dup();
    if (!callable)
        throw logic_error("mutate function must not be NULL");

    return opaque(call_callable(callable, just_kw(kw.dup())));
}

} // namespace sp2
} // namespace python