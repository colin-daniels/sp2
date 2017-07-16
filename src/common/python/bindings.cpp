//
// Created by lampam on 7/11/17.
//

// Python.h must be included prior to any other headers.
#include "Python.h"
#include <common/python/bindings.hpp>

#include <common/python/util.hpp>
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
// helper functions

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
py_scoped_t call_module_function(const char *mod_name, const char *func_name, py_scoped_t args, py_scoped_t kw)
{
    py_scoped_t module = [&] {
        auto name = scope(PyUnicode_DecodeFSDefault(mod_name));
        throw_on_py_err();

        return scope(PyImport_Import(name.raw()));
    }();

    if (!module) {
        throw_on_py_err();
        throw runtime_error("Error loading python module");
    }

    auto func = scope(PyObject_GetAttrString(module.raw(), func_name));
    if (!(func && PyCallable_Check(func.raw()))) {
        throw_on_py_err();
        throw runtime_error("Error accessing python function");
    }

    py_scoped_t retval = scope(PyObject_Call(func.raw(), args.raw(), kw.raw()));
    if (!retval) {
        throw_on_py_err();
        throw runtime_error("Error during Python function call");
    }

    return std::move(retval);
}

// Call a function in a module with *args but no **kw
py_scoped_t call_module_function_args(const char *mod_name, const char *func_name, py_scoped_t args)
{
    auto kw = scope(PyDict_New());
    return move(call_module_function(mod_name, func_name, move(args), move(kw)));
}

// Call a function in a module with **kw but no *args
py_scoped_t call_module_function_kw(const char *mod_name, const char *func_name, py_scoped_t kw)
{
    auto args = scope(PyTuple_New(0));
    return move(call_module_function(mod_name, func_name, move(args), move(kw)));
}

py_scoped_t numpy_array_from_flat_vec(const vector<double> &v, size_t width)
{
    if (v.size() % width != 0) {
        throw logic_error("flat array not divisible by width");
    }

    // copy data into a brand new array object.
    constexpr int ndim = 2;
    npy_intp dims[ndim] = {npy_intp(v.size() / width), npy_intp(width)};
    auto o2 = PyArray_SimpleNew(ndim, dims, NPY_DOUBLE);
    auto o = scope(o2);

    auto data = (double *)PyArray_DATA((PyArrayObject *)o.raw());
    copy(v.begin(), v.end(), data);

    return move(o);
}

// returns an empty string on success,
// returns an error message if the object does not meet specifications,
//  and throws if any other error occurs.
string validate_2d_array_dims(py_scoped_t o, size_t expect_nrow, size_t expect_ncol)
{
    if (!PyArray_Check(o.raw())) {
        return "object is not an array";
    }

    npy_intp ndim = PyArray_NDIM((PyArrayObject *)o.raw());
    throw_on_py_err("error attempting to read numpy array ndim");
    npy_intp *shape = PyArray_SHAPE((PyArrayObject *)o.raw());
    throw_on_py_err("error attempting to read numpy array shape");

    if (ndim != 2) {
        return "expected a 2D numpy array";
    }

    if (expect_nrow != shape[0] || expect_ncol != shape[1]) {
        return "array shape mismatch:\n"
                       " expected: (" + to_string(expect_nrow) + ", " + to_string(expect_ncol) + ")\n"
                       "   actual: (" + to_string(shape[0]) + ", " + to_string(shape[1]) + ")";
    }

    return "";
}

vector<double> flat_vec_from_numpy_array(py_scoped_t o, size_t expect_nrow, size_t expect_ncol)
{
    // Force the array into a contiguous layout if it isn't.
    int min_depth = 0; // ignore
    int max_depth = 0; // ignore
    auto contiguous = scope(PyArray_ContiguousFromAny(o.raw(), NPY_DOUBLE, min_depth, max_depth));
    throw_on_py_err("error making numpy array contiguous");

    string err = validate_2d_array_dims(contiguous.dup(), expect_nrow, expect_ncol);
    if (!err.empty()) {
        throw runtime_error(err);
    }

    auto data = (double *)PyArray_DATA((PyArrayObject *)contiguous.raw());
    throw_on_py_err("error accessing numpy array data");
    size_t size = PyArray_SIZE((PyArrayObject *)contiguous.raw());
    throw_on_py_err("error accessing numpy array size");

    return vector<double>(data, data + size);
}


vector<double> call_2d_vector_function(const char *mod_name, const char *func_name, vector<double> input, size_t width)
{
    if (input.size() % width != 0) {
        throw logic_error("vector not divisible by width");
    }
    size_t height = input.size() / width;

    auto x = numpy_array_from_flat_vec(input, width);
    auto args = scope(Py_BuildValue("(O)", x.raw()));
    throw_on_py_err("Exception constructing python args tuple.");

    auto retval = call_module_function_args(mod_name, func_name, move(args));

    return flat_vec_from_numpy_array(retval.dup(), height, width);
}

void extend_sys_path(const char *dir) {

    // python literal equivalent:
    //
    // if d in set(sys.path):
    //     sys.path.remove(d)
    // sys.path.insert(0, d)

    auto py_dir = scope(PyUnicode_DecodeFSDefault(dir));
    throw_on_py_err("add_to_sys_path: error decoding filesystem path");

    auto sys_path = scope_dup(PySys_GetObject((char*)"path"));
    throw_on_py_err("add_to_sys_path: error reading sys.path");

    // make a set because I can't find where the C API exposes list.__contains__... >_>
    auto sys_path_set = scope(PySet_New(sys_path.raw()));
    throw_on_py_err("add_to_sys_path: error making set from sys.path");

    // remove an existing entry
    int already_there = PySet_Contains(sys_path_set.raw(), py_dir.raw());
    throw_on_py_err("add_to_sys_path: error checking set membership");
    switch (already_there) {
        case 1:
            PyObject_CallMethod(sys_path.raw(), "remove", "o", py_dir.raw());
            throw_on_py_err("add_to_sys_path: error removing existing item");
            break;
        case 0: break;
        default: throw logic_error("unexpected result from __contains__");
    }

    // prepend
    PyList_Insert(sys_path.raw(), 0, py_dir.raw());
    throw_on_py_err("add_to_sys_path: error inserting item");
}


void extend_sys_path(vector<string> dirs) {
    // reverse order so that the prepended results appear in the requested order
    for (auto it = dirs.rbegin(); it != dirs.rend(); it++) {
        extend_sys_path(it->c_str());
    }
}

void initialize(const char *prog) {
    if (prog) {
        // NOTE: Technically, we are leaking this memory by not setting up any sort of
        //       provision to call 'PyMem_RawFree(program)' after 'Py_FinalizeEx()'.
        //       But this is no big deal since this function should only be called once.
        wchar_t *program = Py_DecodeLocale(prog, NULL);
        if (program) {
            Py_SetProgramName(program);
        } else {
            throw runtime_error("Warning: Could not decode program name for python bindings");
        }
    }

    Py_Initialize();

    initialize_numpy();
    throw_on_py_err("error initializing numpy");

    initialize_fake_modules();
}

int finalize() {
    finalize_fake_modules();
    return Py_FinalizeEx();
}

} // namespace impl
} // namespace python
} // namespace sp2

/* --------------------------------------------------------------------- */
// Public API

void sp2::python::initialize(const char *prog) {
    return impl::initialize(prog);
}

int sp2::python::finalize() {
    return impl::finalize();
}

void sp2::python::extend_sys_path(vector<string> dirs) {
    return impl::extend_sys_path(dirs);
}

vector<double> sp2::python::call_2d_vector_function(const char *mod_name, const char *func_name, vector<double> input,
        size_t width) {
    return impl::call_2d_vector_function(mod_name, func_name, input, width);
}

