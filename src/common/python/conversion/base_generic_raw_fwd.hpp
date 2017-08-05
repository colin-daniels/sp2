/* ========================================================================== */
#ifndef SP2_PYTHON_CONVERSION_GENERIC_RAW_FWD
#define SP2_PYTHON_CONVERSION_GENERIC_RAW_FWD

#include "common/python/types/py_ref_t_fwd.hpp"

#undef SP2_PYTHON_CONVERSION_GENERIC_RAW_FWD
#else
#include "diagnostic/forward_dependency_cycle"
#endif
/* ========================================================================== */

#include <vector>
#include <tuple>

namespace sp2 {
namespace python {

template<typename T>
bool to_python(const std::vector<T> &vec, py_ref_t &list);

template<typename T>
bool from_python(const py_ref_t &o, std::vector<T> &vec);

template<typename... Ts>
bool from_python(const py_ref_t &py, std::tuple<Ts...> &cs);

template<typename... Ts>
bool to_python(const std::tuple<Ts...> &cs, py_ref_t &py);

} // namespace python
} // namespace sp2
