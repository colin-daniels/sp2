/* ========================================================================== */
#ifndef SP2_PYTHON_CONVERSION_BASE_MONOMORPHIC_FWD_DAG
#define SP2_PYTHON_CONVERSION_BASE_MONOMORPHIC_FWD_DAG

#include "common/python/types/py_ref_t_fwd.hpp"
#include "common/python/types/py_object_t_fwd.hpp"

#undef SP2_PYTHON_CONVERSION_BASE_MONOMORPHIC_FWD_DAG
#else
#include "diagnostic/forward_dependency_cycle"
#endif
/* ========================================================================== */
#ifndef SP2_PYTHON_CONVERSION_BASE_MONOMORPHIC_FWD
#define SP2_PYTHON_CONVERSION_BASE_MONOMORPHIC_FWD

#include "concrete.hpp"
#include "numpy.hpp"
#include <type_traits>

namespace sp2 {
namespace python {

/* -------------------------------------------------------------------------- */
// Conversions from concrete.hpp,
//  available to all consumers.

#ifdef Py_PYTHON_H
template<
    typename T,
    typename = std::enable_if_t<has_concrete_to_python<T>::value>>
bool to_python(const T &c, py_ref_t &py);
#endif // Py_PYTHON_H

template<
    typename T,
    typename = std::enable_if_t<has_concrete_to_python<T>::value>>
bool to_python(const T &c, py_object_t &py);

#ifdef Py_PYTHON_H
template<
    typename T,
    typename = std::enable_if_t<has_concrete_from_python<T>::value>>
bool from_python(const py_ref_t &py, T &c);
#endif // Py_PYTHON_H

template<
    typename T,
    typename = std::enable_if_t<has_concrete_from_python<T>::value>>
bool from_python(const py_object_t &py, T &c);

/* -------------------------------------------------------------------------- */
// Conversions from numpy.hpp,
//  available to all consumers.

#ifdef Py_PYTHON_H
template<typename T>
bool to_python(const as_ndarray_t<T> &c, py_ref_t &py);
#endif // Py_PYTHON_H

template<typename T>
bool to_python(const as_ndarray_t<T> &c, py_object_t &py);

#ifdef Py_PYTHON_H
template<typename T>
bool from_python(const py_ref_t &py, as_ndarray_t<T> &c);
#endif // Py_PYTHON_H

template<typename T>
bool from_python(const py_object_t &py, as_ndarray_t<T> &c);

} // namespace python
} // namespace sp2

#endif