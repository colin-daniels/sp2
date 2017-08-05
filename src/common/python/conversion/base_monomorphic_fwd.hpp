/* ========================================================================== */
#ifndef SP2_PYTHON_CONVERSION_BASE_MONOMORPHIC_FWD
#define SP2_PYTHON_CONVERSION_BASE_MONOMORPHIC_FWD

#include "common/python/types/py_scoped_t_body_fwd.hpp"
#include "common/python/types/py_opaque_t_body_fwd.hpp"

#undef SP2_PYTHON_CONVERSION_BASE_MONOMORPHIC_FWD
#else
#include "diagnostic/forward_dependency_cycle"
#endif
/* ========================================================================== */

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
    typename = std::enable_if_t<has_concrete_to_python<T>::value>
>
bool to_python(const T &c, py_scoped_t &py);
#endif // Py_PYTHON_H

template<
    typename T,
    typename = std::enable_if_t<has_concrete_to_python<T>::value>
>
bool to_python(const T &c, py_opaque_t &py);

#ifdef Py_PYTHON_H
template<
    typename T,
    typename = std::enable_if_t<has_concrete_from_python<T>::value>
>
bool from_python(const py_scoped_t &py, T &c);
#endif // Py_PYTHON_H

template<
    typename T,
    typename = std::enable_if_t<has_concrete_from_python<T>::value>
>
bool from_python(const py_opaque_t &py, T &c);

/* -------------------------------------------------------------------------- */
// Conversions from numpy.hpp,
//  available to all consumers.

#ifdef Py_PYTHON_H
template<typename T>
bool to_python(const as_ndarray_t<T> &c, py_scoped_t &py);
#endif // Py_PYTHON_H

template<typename T>
bool to_python(const as_ndarray_t<T> &c, py_opaque_t &py);

#ifdef Py_PYTHON_H
template<typename T>
bool from_python(const py_scoped_t &py, as_ndarray_t<T> &c);
#endif // Py_PYTHON_H

template<typename T>
bool from_python(const py_opaque_t &py, as_ndarray_t<T> &c);

} // namespace python
} // namespace sp2