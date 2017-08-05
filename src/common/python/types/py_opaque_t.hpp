#ifndef SP2_PY_OPAQUE_T_HPP
#define SP2_PY_OPAQUE_T_HPP

#include "py_opaque_t_body.hpp"
#include "common/python/conversion.hpp"

#include <string>

namespace sp2 {
namespace python {

template<typename T>
py_opaque_t py_opaque_t::from(const T &c, const char *msg)
{
    py_opaque_t py;
    if (!to_python(c, py))
        throw std::runtime_error(msg);
    return py;
}

template<
    typename T,
    typename = std::enable_if_t<std::is_default_constructible<T>::value>
>
T py_opaque_t::parse_as(const char *msg) const
{
    auto c = T();
    if (!from_python(*this, c))
        throw std::runtime_error(msg);
    return c;
}

} // namespace python
} // namespace sp2
#endif // header guard