#ifndef SP2_PYTHON_INTERNALS_HPP
#define SP2_PYTHON_INTERNALS_HPP

#include "diagnostic/expect_python_headers"

// General include file for utilities that assist in direct interaction
// with the CPython API.
//
// CPP files in common/python should include this unless they know better.
// Also, their very first #include must be <Python.h>.

// You get:

// * py_ref_t
//   An RAII wrapper type around Python references that simplifies reasoning
//   about refcounts and exception safety.
//
// * from_python, to_python
//   Conversion functions which will help eliminate some API calls
//
// * print_on_py_err, throw_on_py_err
//   These help partially automate the process of checking the PyErr API.

#include "types/py_ref_t.hpp"
#include "conversion.hpp"
#include "error.hpp"

#endif // header guard