#include <Python.h> // Must be first include

#include "error.hpp"

#include <stdexcept>

namespace sp2 {
namespace python {

void throw_on_py_err(const char *msg)
{
    if (PyErr_Occurred())
    {
        PyErr_Print();
        throw std::runtime_error(msg);
    }
}

void throw_on_py_err()
{
    throw_on_py_err("An exception was thrown in Python.");
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

} // namespace python
} // namespace sp2