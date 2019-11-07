#ifndef SP2_PYTHON_UTILITY_HPP
#define SP2_PYTHON_UTILITY_HPP

#include "common/python/types/py_object_t_fwd.hpp"

namespace sp2 {
namespace python {

/// Conflict-resolution strategies for merge_dictionaries.
enum class merge_strategy : int
{
    /// Resolve key conflicts by taking the first dict's value.
        USE_FIRST = 0,
    /// Resolve key conflicts by taking the second dict's value.
        USE_SECOND = 1,
    /// Don't resolve key conflicts; throw a runtime_exception
        ERROR = 2
};

/// Perform a union-like operation on two python dictionaries that produces a
/// dict with all of their (key, value) pairs.
py_object_t merge_dictionaries(
    const py_object_t &a, const py_object_t &b,
    merge_strategy strategy = merge_strategy::USE_SECOND);


} // namespace sp2
} // namespace python


#endif // SP2_PYTHON_UTILITY_HPP
