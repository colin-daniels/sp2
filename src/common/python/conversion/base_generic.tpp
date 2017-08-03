/* ========================================================================== */
#ifndef SP2_PYTHON_CONVERSION_BASE_GENERIC_TPP
#define SP2_PYTHON_CONVERSION_BASE_GENERIC_TPP

#include "common/python/py_opaque_t.tpp"

#undef SP2_PYTHON_CONVERSION_BASE_GENERIC_TPP
#else
#include "diagnostic/forward_dependency_cycle"
#endif // SP2_PYTHON_CONVERSION_BASE_GENERIC_TPP
/* ========================================================================== */

#include <vector>
#include <tuple>

namespace sp2 {
namespace python {

template<typename T>
bool to_python(const std::vector<T> &vec, py_opaque_t &py_list);

template<typename T>
bool from_python(const py_opaque_t &o, std::vector<T> &vec);

template<typename... Ts>
bool from_python(const py_opaque_t &py, std::tuple<Ts...> &cs);

template<typename... Ts>
bool to_python(const std::tuple<Ts...> &cs, py_opaque_t &py);

} // namespace python
} // namespace sp2
