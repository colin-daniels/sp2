#ifndef SP2_PYTHON_PUBLIC_HPP
#define SP2_PYTHON_PUBLIC_HPP

// General include file for the outward-facing API of the python directory.
//
// This is the recommended file to import for CPP files outside of
// common/python.  Headers can also import this, or choose to include smaller
// pieces of it at their discretion.
//
// Some more lower-level functionality can be found
//  in "common/python/internals.hpp".

// For including this, you get:
//
// * environment
//   for initialization and cleanup of the interpreter
//
// * py_object_t
//   a wrapper type around python references with 'getattr' and 'call' member
//   functions (enabling arbitrary interaction with them)
//
// * py_import("name"), py_from(x), py_tuple(...)
//   for creating objects
//
// * the parse<T>() member function
//   for reading values back into C++ data types
//
// * py_error
//   a C++ exception that represents a Python exception,
//   complete with backtrace
//
// * sp2::python::fake_modules
//   Direct access to the singleton instances of python modules which are
//   embedded into the C++ source code.

#include "common/python/types/py_object_t.hpp"
#include "common/python/environment.hpp"
#include "common/python/error.hpp"
#include "common/python/modules/fake_modules.hpp"

#endif // header guard