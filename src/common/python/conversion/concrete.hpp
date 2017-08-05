#ifndef SP2_PYTHON_CONVERSION_CONCRETE_HPP
#define SP2_PYTHON_CONVERSION_CONCRETE_HPP

#include "common/python/types/py_ref_t_fwd.hpp"
#include "common/python/types/py_object_t_fwd.hpp"
#include "common/minimize/metropolis_enums.hpp"

#include <string>
#include <type_traits>

namespace sp2 {
namespace python {

/* -------------------------------------------------------------------------- */
// Macros to iterate over types implemented through `{to,from}_python_concrete`.
// (note: as_ndarray is also supported, but that has its own iteration macro)

// List of types implemented in 'primitive.cpp'
#define SP2_FOR_EACH_CONCRETE_PRIMITIVE(mac) \
mac(long); \
mac(long long); \
mac(unsigned long); \
mac(unsigned long long); \
mac(double); \
mac(bool); \
mac(std::string); \
mac(std::nullptr_t); \
mac(sp2::python::py_ref_t); \
mac(sp2::python::py_object_t)

// Lists of types implemented elsewhere; possibly only in one direction.
#define SP2_FOR_EACH_CONCRETE_CUSTOM_FROM(mac) \
mac(sp2::structural_mutation_type); \
mac(sp2::structural_mutation_t)

#define SP2_FOR_EACH_CONCRETE_CUSTOM_TO(mac)

// List of all types with concrete implementations.
#define SP2_FOR_EACH_CONCRETE_FROM(mac) \
SP2_FOR_EACH_CONCRETE_PRIMITIVE(mac); \
SP2_FOR_EACH_CONCRETE_CUSTOM_FROM(mac)

#define SP2_FOR_EACH_CONCRETE_TO(mac) \
SP2_FOR_EACH_CONCRETE_PRIMITIVE(mac); \
SP2_FOR_EACH_CONCRETE_CUSTOM_TO(mac)

/* -------------------------------------------------------------------------- */
// actually generate the declarations

#define MAC(T) \
    bool to_python_concrete(const T &c, py_ref_t &py)
SP2_FOR_EACH_CONCRETE_TO( MAC );
#undef MAC

#define MAC(T) \
    bool from_python_concrete(const py_ref_t &py, T &c)
SP2_FOR_EACH_CONCRETE_FROM( MAC );
#undef MAC

/* -------------------------------------------------------------------------- */
// type trait for having a concrete implementation

template<typename T>
struct has_concrete_to_python : public std::false_type {};
template<typename T>
struct has_concrete_from_python : public std::false_type {};

// FIXME were my SFINAE-foo a little better, I could probably find another way
//       to do this beyond just generating a ton of specializations.

#define MAC(T) \
    template<> \
    struct has_concrete_from_python<T> : public std::true_type {}

SP2_FOR_EACH_CONCRETE_FROM(MAC);
#undef MAC

#define MAC(T) \
    template<> \
    struct has_concrete_to_python<T> : public std::true_type {}

SP2_FOR_EACH_CONCRETE_TO(MAC);
#undef MAC

} // namespace python
} // namespace sp2

#endif // SP2_PYTHON_CONVERSION_CONCRETE_HPP
