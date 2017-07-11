//
// Created by lampam on 7/11/17.
//

#include "bindings.h"

using namespace std;
using namespace sp2::python;

sp2::python::py_scoped_t::py_scoped_t(PyObject *o)
        : obj(o)
{ }

sp2::python::py_scoped_t sp2::python::scope_dup(PyObject *o) {
    Py_XINCREF(o);
    return py_scoped_t(o);
}

sp2::python::py_scoped_t sp2::python::scope(PyObject *o) {
    return py_scoped_t(o);
}

void ::sp2::python::throw_on_py_err(const char *msg) {
    if (PyErr_Occurred())
        PyErr_Print();
    throw runtime_error(msg);
}

void ::sp2::python::throw_on_py_err() {
    throw_on_py_err("An exception was thrown in Python.");
}

void ::sp2::python::initialize(const char *prog) {
    if (prog) {
        wchar_t *program = Py_DecodeLocale(prog, NULL);
        if (program) {
            Py_SetProgramName(program);
        } else {
            throw runtime_error("Warning: Could not decode program name for python bindings");
        }
    }

    Py_Initialize();
}

int ::sp2::python::finalize() {
    return Py_FinalizeEx();
}

/*
std::vector<double> sp2::python::flat_2d_from_py_sequence(size_t ncol, PyObject *o )
{
    if (!PySequence_Check(o)) {
        throw runtime_error("Python value was not a sequence");
    }

    Py_ssize_t len = PySequence_Length(o);
    if (len == -1) {
        throw runtime_error("Could not get python sequence length");
    }

    vector<double> positions(ncol * len);
    for (int i=0; i < len; i++) {
        auto item = scope(PySequence_GetItem(o, i));

        auto subsequence = PySequence_Fast(item.raw(), "noniterable subsequence");
        throw_on_py_err();
        if (!subsequence) {
            throw runtime_error("Unknown error constructing sequence");
        }

        if (ncol != PySequence_Length(subsequence.raw())) {
            throw runtime_error("Expected subsequences of fixed length");
        }

        for (int k=0; k < ncol; k++) {
            // borrowed reference; don't decref
            PyObject *pItem = PySequence_Fast_GET_ITEM(subsequence.raw(), k);

            subsequence[ncol*i + k] = PyFloat_AsDouble(pItem);
            throw_on_py_err();
        }
    }

    return positions;
}
 */

std::vector<double> sp2::python::vec_from_py_sequence(PyObject *o)
{
    if (!PySequence_Check(o)) {
        throw runtime_error("Python value was not a sequence");
    }

    Py_ssize_t len = PySequence_Length(o);
    if (len == -1) {
        throw_on_py_err();
        throw runtime_error("Could not get python sequence length");
    }

    vector<double> out(len);
    for (int i=0; i < len; i++) {
        auto item = scope(PySequence_GetItem(o, i));
        throw_on_py_err();

        out[i] = PyFloat_AsDouble(item.raw());
        throw_on_py_err();
    }

    return out;
}

sp2::python::py_scoped_t sp2::python::call_module_function(const char *mod_name, const char *func_name,
                                                            sp2::python::py_scoped_t args, sp2::python::py_scoped_t kw) {
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
        throw runtime_error("Error getting Python retval");
    }

    return std::move(retval);
}

sp2::python::py_scoped_t
sp2::python::call_module_function_kw(const char *mod_name, const char *func_name,
                                     py_scoped_t kw) {
    auto args = scope(PyTuple_New(0));
    return std::move(call_module_function(mod_name, func_name, move(args), move(kw)));
}

std::vector<double> sp2::python::call_vector_function(const char *mod_name, const char *func_name, std::vector<double> input) {
    auto kw = scope(PyDict_New());
    if (!kw) {
        throw_on_py_err();
        throw runtime_error("Exception constructing python dictionary.");
    }

    {
        auto x = py_list_from_vec(input);
        if (PyDict_SetItemString(kw.raw(), "x", x.raw())) {
            throw_on_py_err();
            throw runtime_error("Exception setting python dict item.");
        }
    }

    auto retval = call_module_function_kw(mod_name, func_name, move(kw));

    return vec_from_py_sequence(retval.raw());
}

py_scoped_t sp2::python::py_list_from_vec(std::vector<double> v)
{
    auto list = scope(PyList_New(v.size()));
    if (!list) {
        throw new runtime_error("error constructing python list");
    }
    for (int i=0; i < v.size(); i++) {
        /* FIXME
         * if this throws then 'list' will get its refcount decremented to zero while
         *  it still holds some NULL elements. Is the list destructor robust to this?
         */
        auto value = scope(PyFloat_FromDouble(v[i]));
        throw_on_py_err();

        PyList_SET_ITEM(list.raw(), i, value.steal());
        throw_on_py_err();
    }
    return list;
}

// move constructor
sp2::python::py_scoped_t::py_scoped_t(sp2::python::py_scoped_t &&other)
        : obj(other.steal())
{ }

// move assignment operator
sp2::python::py_scoped_t &sp2::python::py_scoped_t::operator=(sp2::python::py_scoped_t &&other) {
    if (this != &other) {
        obj = other.steal();
    }
    return *this;
}

sp2::python::py_scoped_t::~py_scoped_t() {
    Py_XDECREF(obj);
}

sp2::python::py_scoped_t sp2::python::py_scoped_t::dup() {
    Py_XINCREF(obj);
    return py_scoped_t(obj);
}

PyObject *sp2::python::py_scoped_t::raw() const {
    return obj;
}

PyObject *sp2::python::py_scoped_t::steal() {
    auto tmp = obj;
    obj = NULL;
    return tmp;
}

sp2::python::py_scoped_t::operator bool() const {
    return (bool)obj;
}
