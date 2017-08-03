#ifndef SP2_PYTHON_CONVERSION_HPP
#define SP2_PYTHON_CONVERSION_HPP

// These files provide {to,from}_python.
// see common/python/conversion/README.md for details.

#include "common/python/conversion/base_monomorphic.hpp"
#include "common/python/conversion/base_generic.hpp"

#ifdef Py_PYTHON_H
// a special treat for those who make the sacrifice of importing Python.h...
#include "common/python/conversion/base_generic_raw.hpp"
#endif // Py_PYTHON_H

#endif // SP2_PYTHON_CONVERSION_HPP
