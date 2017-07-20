#include "Python.h" // Must be first include

#include <common/python/util.hpp>

#include <common/util/templates.hpp>

#include <stdexcept>

using namespace std;

namespace sp2 {
namespace python {

py_scoped_t::py_scoped_t(PyObject *o)
        : obj(o)
{ }

py_scoped_t::operator bool() const {
    return (bool)obj;
}

py_scoped_t::py_scoped_t(py_scoped_t &&other)
        : obj(other.steal())
{ }

py_scoped_t &py_scoped_t::operator=(py_scoped_t &&other) {
    if (this != &other) {
        if (obj) {
            // we could implicitly destroy the existing reference, but with no
            // clear use case, the conservative choice is to require explicit
            // destruction
            throw logic_error("attempted to overwrite occupied py_scoped_t");
        }
        obj = other.steal();
    }
    return *this;
}

py_scoped_t::~py_scoped_t() {
    destroy();
}

py_scoped_t py_scoped_t::dup() {
    Py_XINCREF(obj);
    return py_scoped_t(obj);
}

PyObject *py_scoped_t::raw() {
    return obj;
}

PyObject *py_scoped_t::steal() {
    auto tmp = obj;
    obj = NULL;
    return tmp;
}

void py_scoped_t::destroy() {
    Py_XDECREF(obj);
    obj = NULL;
}

py_scoped_t scope(PyObject *o) {
    return py_scoped_t(o);
}

py_scoped_t scope_dup(PyObject *o) {
    Py_XINCREF(o);
    return py_scoped_t(o);
}

// --------------------------------

void throw_on_py_err(const char *msg) {
    if (PyErr_Occurred()) {
        PyErr_Print();
        throw runtime_error(msg);
    }
}

void throw_on_py_err() {
    throw_on_py_err("An exception was thrown in Python.");
}

// --------------------------------

// Implementation of repr() and str().
// F is a function(PyObject *) -> PyObject * returning a new reference to a unicode 'str' object
template <typename F>
std::string str_impl(py_scoped_t &o, F stringify) {

    auto py_str = scope(stringify(o.raw()));
    throw_on_py_err("repr: error stringifying");

    // NOTE: we are not responsible for deallocating this;
    //       it's lifetime is bound to the py_str python object.
    char *str = PyUnicode_AsUTF8(py_str.raw());
    throw_on_py_err("repr: error encoding as utf8");
    auto guard = sp2::scope_guard([&] { PyMem_Free(str); });

    return string(str);
}

string str(py_scoped_t &o) {
    return str_impl(o, [](auto x) { return PyObject_Str(x); });
}

string repr(py_scoped_t &o) {
    return str_impl(o, [](auto x) { return PyObject_Repr(x); });
}

bool print_on_py_err() {
    if (PyErr_Occurred()) {
        PyErr_Print(); // Note: clears the error.
        return true;
    }
    return false;
}

} // namespace python
} // namespace sp2
