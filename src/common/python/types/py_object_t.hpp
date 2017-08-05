#ifndef SP2_PY_OBJECT_T_HPP
#define SP2_PY_OBJECT_T_HPP

#include "py_object_t_fwd.hpp"

#include "common/python/conversion.hpp"

#include <string>

#include "py_object_t_body.hpp"

namespace sp2 {
namespace python {

/// Value-returning wrapper around to_python.
///
/// Errors are communicated by std::runtime_error.
template<typename T>
py_object_t py_from(const T &c, const char *msg)
{
    py_object_t py;
    if (!to_python(c, py))
        throw std::runtime_error(msg);
    return py;
}

template<typename T>
py_object_t py_from(const T &c)
{ return py_from(c, "an error occurred converting data into python objects"); }

template<typename... Ts>
py_object_t py_tuple(Ts... ts)
{ return detail::py_tuple_noconv({py_from(ts)...}); }

template<
    typename T,
    typename = std::enable_if_t<std::is_default_constructible<T>::value>
>
T py_object_t::parse_as(const char *msg) const
{
    auto c = T();
    if (!from_python(*this, c))
        throw std::runtime_error(msg);
    return c;
}

} // namespace python
} // namespace sp2
#endif // header guard