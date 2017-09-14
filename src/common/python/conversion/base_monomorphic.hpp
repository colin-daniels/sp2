#ifndef SP2_PYTHON_CONVERSION_BASE_MONOMORPHIC_HPP
#define SP2_PYTHON_CONVERSION_BASE_MONOMORPHIC_HPP

#include "base_monomorphic_fwd.hpp"

#include "concrete.hpp"
#include "numpy.hpp"
#include "common/python/types/py_object_t_body.hpp"

#include <type_traits>

// base_monomorphic.hpp:
//
//  Provides '{to,from}_python' that forward to functions that have an explicit
//  set of overloads compiled into the binary.

namespace sp2 {
namespace python {

/* -------------------------------------------------------------------------- */
// Conversions from concrete.hpp,
//  available to all consumers.

#ifdef Py_PYTHON_H
template<typename T, typename>
bool to_python(const T &c, py_ref_t &py)
{ return to_python_concrete(c, py); }
#endif // Py_PYTHON_H

template<typename T, typename>
bool to_python(const T &c, py_object_t &py)
{ return to_python_concrete(c, py.inner()); }

#ifdef Py_PYTHON_H
template<typename T, typename>
bool from_python(const py_ref_t &py, T &c)
{ return from_python_concrete(py, c); }
#endif // Py_PYTHON_H

template<typename T, typename>
bool from_python(const py_object_t &py, T &c)
{ return from_python_concrete(py.inner(), c); }

/* -------------------------------------------------------------------------- */
// Conversions from numpy.hpp,
//  available to all consumers.

#ifdef Py_PYTHON_H
template<typename T>
bool to_python(const as_ndarray_t<T> &c, py_ref_t &py)
{ return to_python_by_ndarray(c, py); }
#endif // Py_PYTHON_H

template<typename T>
bool to_python(const as_ndarray_t<T> &c, py_object_t &py)
{ return to_python_by_ndarray(c, py.inner()); }

#ifdef Py_PYTHON_H
template<typename T>
bool from_python(const py_ref_t &py, as_ndarray_t<T> &c)
{ return from_python_by_ndarray(py, c); }
#endif // Py_PYTHON_H

template<typename T>
bool from_python(const py_object_t &py, as_ndarray_t<T> &c)
{ return from_python_by_ndarray(py.inner(), c); }

} // namespace python
} // namespace sp2

#endif // SP2_PYTHON_CONVERSION_BASE_MONOMORPHIC_HPP
