#ifndef SP2_PYTHON_CONVERSION_TPP
#define SP2_PYTHON_CONVERSION_TPP

// These files provide {to,from}_python.
// see common/python/conversion/README.md for details.

#include "common/python/conversion/base_monomorphic.tpp"
#include "common/python/conversion/base_generic.tpp"

#ifdef Py_PYTHON_H
// a special treat for those who make the sacrifice of importing Python.h...
#include "common/python/conversion/base_generic_raw.tpp"
#endif // Py_PYTHON_H

#endif // SP2_PYTHON_CONVERSION_TPP
