#ifndef SP2_UTIL_HPP
#define SP2_UTIL_HPP

#ifndef Py_PYTHON_H
#error This module is designed for internal use by source files in common/python, \
and is not safe to use unless Python.h is already included before all standard \
library headers.
#endif // Py_PYTHON_H

#include <string>
#include <utility>
#include <stdexcept>

namespace sp2 {
namespace python {

// A scoped reference to a python object that uses RAII to handle Py_DECREF.
// This makes it somewhat easier to reason about exception safety,
// though it is not a panacea.
//
// One should still be careful to consider destruction order (since a decref can
//  potentially invoke arbitrary python code), and read the Python API docs
//  carefully to understand when references are duplicated, borrowed, and stolen.
//
// The default PyObject * constructor does NOT perform an incref, since the
// majority of Python API functions return a freshly-incremented reference.
// For those rare functions that return borrowed references, you should
// use the explicit 'scope_dup' constructor instead.
//
// The contained object may be NULL.
class py_scoped_t
{
    PyObject *obj = nullptr;

public:

    // null constructor
    py_scoped_t() {};

    // PyObject constructor
    explicit py_scoped_t(PyObject * o);

    explicit operator bool() const;

    // No copying.  Use dup() to make the incref explicit.
    py_scoped_t& operator=(const py_scoped_t& other) = delete;
    py_scoped_t(const py_scoped_t& other) = delete;

    // move constructor
    py_scoped_t(py_scoped_t&& other);

    // move assignment operator
    py_scoped_t& operator=(py_scoped_t&& other);

    ~py_scoped_t();

    // Increment the refcount and return a new scoped reference.
    py_scoped_t dup();

    // Borrow the reference without touching the refcount.
    //
    // This is the appropriate method for interfacing with most Python APIs.
    PyObject * raw();

    // Leak the reference, preventing the DECREF that would otherwise occur at scope exit.
    // The scoped reference will become NULL.
    //
    // Necessary for working with Python API functions that steal references,
    // such as PyTuple_SetItem.
    PyObject * steal();

    // Explicit destructor.
    //
    // Destroy the reference early, decrementing the refcount and nulling out the pointer
    //  so that nothing happens at scope exit.
    //
    // This can be used to explicitly control the destruction order in places where
    // the natural order of destruction would not be safe.
    void destroy();
};

// explicit constructor from a new ref
py_scoped_t scope(PyObject * o);

// explicit constructor from a borrowed ref, which makes a new reference
py_scoped_t scope_dup(PyObject * o);

// --------------------------------

// Checks for python exceptions via the PyErr API and turns them into C++ exceptions.
void throw_on_py_err(const char *msg);

void throw_on_py_err();

// --------------------------------

// get an object's repr(), mostly for debug purposes
std::wstring repr(py_scoped_t o);

// get an object's str(), mostly for debug purposes
std::wstring str(py_scoped_t o);

} // namespace python
} // namespace sp2

#endif //SP2_UTIL_HPP
