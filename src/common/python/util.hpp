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

// The designated owner of a Python object reference.
//
// This scoped reference to a python object that uses RAII to handle Py_DECREF.
// This makes it somewhat easier to reason about exception safety,
// though it is not a panacea.
//
// One should still be careful to consider destruction order (since a decref can
// potentially invoke arbitrary python code), and read the Python API docs
// carefully to understand when references are duplicated, borrowed, and stolen.
//
// The default PyObject * constructor does NOT perform an incref, since the
// majority of Python API functions return a freshly-incremented reference.
// For those rare functions that return borrowed references, you should
// use the explicit 'scope_dup' constructor instead.
//
// 'const' guarantees for this type are particularly weak. In general, a
// 'const py_scoped_t' won't allow you to modify its pointer to point somewhere
// else, but you are free to modify its referent in any other way.
// This admission is made because a 'const PyObject*' is nearly useless.
//
// The contained object may be NULL.
class py_scoped_t
{
    mutable PyObject *obj = nullptr;

public:

    // null constructor
    py_scoped_t() {};

    // PyObject constructor
    explicit py_scoped_t(PyObject *o);

    explicit operator bool() const;

    // No copying.  Use dup() to make the incref explicit.
    py_scoped_t &operator=(const py_scoped_t &other) = delete;

    py_scoped_t(const py_scoped_t &other) = delete;

    // move constructor
    py_scoped_t(py_scoped_t &&other);

    // move assignment operator
    py_scoped_t& operator=(py_scoped_t &&other);

    ~py_scoped_t();

    // Increment the refcount and return a new scoped reference.
    py_scoped_t dup() const;

    // Borrow the reference without touching the refcount.
    //
    // This is the appropriate method for interfacing with most Python APIs.
    //
    // For reasons discussed earlier, the returned pointer is mutable.
    // Just... be reasonable, okay?
    PyObject *raw() const;

    // Leak the reference, preventing the DECREF that would otherwise occur at scope exit.
    // The scoped reference will become NULL.
    //
    // Necessary for working with Python API functions that steal references,
    // such as PyTuple_SetItem.
    PyObject *steal();

    // Postfix move() for convenience, and for easier interchange with dup().
    //
    // Note a small semantic difference in that, due to the by-value return
    // type, this is guaranteed to clear the receiver.
    py_scoped_t move();

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
py_scoped_t scope(PyObject *o);

// explicit constructor from a borrowed ref, which makes a new reference
py_scoped_t scope_dup(PyObject *o);

// --------------------------------

// Checks for python exceptions via the PyErr API and turns them into C++ exceptions.
void throw_on_py_err(const char *msg);

void throw_on_py_err();

// Checks for a python exception and prints it to stderr, clearing the error.
//
// Returns true if an error was found.
bool print_on_py_err();

// --------------------------------

/// get an object's repr() in utf8, mostly for debug purposes.
std::string repr(py_scoped_t &o);

/// get an object's str() in utf8, mostly for debug purposes.
std::string str(py_scoped_t &o);

// --------------------------------

/// Access an attribute of a python object.
///
/// Equivalent to 'getattr(obj, name)'.
/// Throws an exception if the attribute does not exist.
py_scoped_t getattr(py_scoped_t &o, const char *attr);

/// Access an attribute of a python object, or a default value.
///
/// Equivalent to 'getattr(obj, name, def)'.
/// Throws an exception if the attribute does not exist.
py_scoped_t getattr(py_scoped_t &o, const char *attr, py_scoped_t &def);
py_scoped_t getattr(py_scoped_t &o, const char *attr, py_scoped_t &&def);




} // namespace python
} // namespace sp2

#endif //SP2_UTIL_HPP
