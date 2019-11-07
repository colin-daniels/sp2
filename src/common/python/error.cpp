#include <Python.h> // Must be first include

#include "error.hpp"

#include "common/python/types/py_ref_t.hpp"

#include <stdexcept>

namespace sp2 {
namespace python {

using namespace std;

void throw_on_py_err(const char *msg)
{
    if (PyErr_Occurred())
        throw *grab_py_err(msg);
}

bool print_on_py_err()
{
    if (PyErr_Occurred())
    {
        PyErr_Print(); // Note: clears the error.
        return true;
    }
    return false;
}

std::unique_ptr<py_error> grab_py_err(const char* msg)
{
    PyObject *raw_type, *raw_value, *raw_traceback;
    PyErr_Fetch(&raw_type, &raw_value, &raw_traceback);

    // "If the error indicator is not set, set all three variables to NULL.
    //  If it is set, it will be cleared and you own a reference to each
    //  object retrieved.
    //  The value and traceback object may be NULL even when
    //  the type object is not."
    // - https://docs.python.org/3/c-api/exceptions.html#c.PyErr_Fetch
    if (!raw_type)
        return {};

    py_object_t type = opaque(scope(raw_type));
    py_object_t value = opaque(scope(raw_value));
    py_object_t traceback = opaque(scope(raw_traceback));
    return make_unique<py_error>(py_exc_info{type, value, traceback}, string(msg));
}

void restore_py_err(const py_error &err)
{
    // PyErr_Restore steals references.
    PyObject* raw_type = err.exc_type().inner().dup().steal();
    PyObject* raw_value = err.exc_value().inner().dup().steal();
    PyObject* raw_traceback = err.exc_traceback().inner().dup().steal();
    PyErr_Restore(raw_type, raw_value, raw_traceback);
}

const py_error& py_error::print_to_stderr() const
{
    unique_ptr<py_error> temp = grab_py_err();

    restore_py_err(*this);
    PyErr_Print();

    if (temp)
        restore_py_err(*temp);

    return *this;
}

void py_error::rethrow() const
{
    throw static_cast<std::runtime_error>(*this);
}
} // namespace python
} // namespace sp2