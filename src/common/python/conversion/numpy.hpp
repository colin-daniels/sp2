#ifndef SP2_PYTHON_CONVERSION_NUMPY_HPP
#define SP2_PYTHON_CONVERSION_NUMPY_HPP

#include "common/python/bindings.hpp"
#include "common/python/types/as_ndarray.hpp"
#include "common/python/types/py_ref_t_fwd.hpp"

#include <cstdint>
#include <vector>
#include <string>

namespace sp2 {
namespace python {

/* -------------------------------------------------------------------------- */
// Conversions of as_ndarray_t <-> ndarray
//
// These only support a fixed set of primitive dtypes.
// Technically, arrays of arbitrary objects do exist,
// but they are terrible footguns.

#define SP2_FOR_EACH_NUMPY_DTYPE(mac) \
mac(bool); \
mac(std::int8_t); \
mac(std::int16_t); \
mac(std::int32_t); \
mac(std::int64_t); \
mac(std::uint8_t); \
mac(std::uint16_t); \
mac(std::uint32_t); \
mac(std::uint64_t); \
mac(float); \
mac(double)

//--------
// Forward declare the generic functions used for implementation.
/// Implementation detail.  Use `to_python` instead.
template<typename T>
bool to_python_by_ndarray(const as_ndarray_t<T> &c, py_ref_t &py);

/// Implementation detail.  Use `from_python` instead.
template<typename T>
bool from_python_by_ndarray(const py_ref_t &py, as_ndarray_t<T> &c);

//--------
// Declare all instantiations, EXCEPT where they are actually defined.
#ifndef SP2_NDARRAY_INSTANTIATIONS_VISIBLE

#define MAC(T) \
extern template bool to_python_by_ndarray<T>(const as_ndarray_t<T> &c, py_ref_t &py); \
extern template bool from_python_by_ndarray<T>(const py_ref_t &py, as_ndarray_t<T> &c)

SP2_FOR_EACH_NUMPY_DTYPE(MAC);
#undef MAC

#endif // !SP2_NDARRAY_INSTANTIATIONS_VISIBLE

} // namespace python
} // namespace sp2

#endif // SP2_PYTHON_CONVERSION_NUMPY_HPP
