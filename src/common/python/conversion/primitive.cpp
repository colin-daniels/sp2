#include <Python.h> // must be first include

#include "concrete.hpp"

#include "common/python/internals.hpp"

using namespace std;

namespace sp2 {
namespace python {

/* --------------------------------------------------------------------- */

void refuse_to_evict(py_ref_t &p)
{
    if (p)
    {
        throw std::logic_error(
            "attempted to serialize into already-occupied py_ref_t");
    }
}

// basic implementation given a function with these properties:
//  - type (Args...) -> PyObject*
//  - Returns a new reference on success
//  - Returns NULL and sets the python error state on failure
template<typename F, typename... Args>
bool to_python_simple_conv(py_ref_t &py, F func, Args... args)
{
    refuse_to_evict(py);

    py = scope(func(args...));
    return !print_on_py_err();
}

bool to_python_concrete(const nullptr_t &c, py_ref_t &py)
{
    refuse_to_evict(py);

    py = scope_dup(Py_None);
    return true;
}

bool to_python_concrete(const long &c, py_ref_t &py)
{
    return to_python_simple_conv(py, PyLong_FromLong, c);
}

bool to_python_concrete(const unsigned long &c, py_ref_t &py)
{
    return to_python_simple_conv(py, PyLong_FromUnsignedLong, c);
}

bool to_python_concrete(const long long &c, py_ref_t &py)
{
    return to_python_simple_conv(py, PyLong_FromLongLong, c);
}

bool to_python_concrete(const unsigned long long &c, py_ref_t &py)
{
    return to_python_simple_conv(py, PyLong_FromUnsignedLongLong, c);
}

bool to_python_concrete(const double &c, py_ref_t &py)
{
    return to_python_simple_conv(py, PyFloat_FromDouble, c);
}

bool to_python_concrete(const bool &c, py_ref_t &py)
{
    refuse_to_evict(py);

    if (c)
        py = scope_dup(Py_True);
    else
        py = scope_dup(Py_False);

    return true;
}

bool to_python_concrete(const std::string &utf8, py_ref_t &py)
{
    return to_python_simple_conv(py,
        PyUnicode_DecodeUTF8, utf8.c_str(), utf8.size(), "strict");
}

bool to_python_concrete(const py_ref_t &c, py_ref_t &py)
{
    refuse_to_evict(py);

    py = c.dup();
    return true;
}

bool to_python_concrete(const py_object_t &c, py_ref_t &py)
{
    refuse_to_evict(py);

    py = c.inner();
    return true;
}

/* -------------------------------------------------------------------------- */

bool from_python_concrete(const py_ref_t &py, long &c)
{
    c = PyLong_AsLong(py.raw());
    return !print_on_py_err();
}

bool from_python_concrete(const py_ref_t &py, unsigned long &c)
{
    c = PyLong_AsUnsignedLong(py.raw());
    return !print_on_py_err();
}

bool from_python_concrete(const py_ref_t &py, long long &c)
{
    c = PyLong_AsLongLong(py.raw());
    return !print_on_py_err();
}

bool from_python_concrete(const py_ref_t &py, unsigned long long &c)
{
    c = PyLong_AsUnsignedLongLong(py.raw());
    return !print_on_py_err();
}

bool from_python_concrete(const py_ref_t &py, double &c)
{
    c = PyFloat_AsDouble(py.raw());
    return !print_on_py_err();
}

bool from_python_concrete(const py_ref_t &py, bool &c)
{
    if (py.raw() == Py_True)
    {
        c = true;
        return true;

    }
    else if (py.raw() == Py_False)
    {
        c = false;
        return true;
    }

    return false;
}

bool from_python_concrete(const py_ref_t &py, nullptr_t &c)
{
    c = nullptr; // *shrug*
    return py.raw() == Py_None;
}

bool from_python_concrete(const py_ref_t &py, std::string &c)
{
    c.clear();

    // From the Python API docs:
    //     This caches the UTF-8 representation of the string in the
    //     Unicode object, and subsequent calls will return a pointer
    //     to the same buffer.
    //     The caller is not responsible for deallocating the buffer.
    char *buf = PyUnicode_AsUTF8(py.raw());
    if (!buf && print_on_py_err())
    {
        return false;
    }
    c = std::string(buf);
    return true;
}

bool from_python_concrete(const py_ref_t &py, py_ref_t &c)
{
    refuse_to_evict(c);

    c = py.dup();
    return true;
}

bool from_python_concrete(const py_ref_t &py, py_object_t &c)
{
    refuse_to_evict(c.inner());

    c.inner() = py.dup();
    return true;
}

} // namespace python
} // namespace sp2
