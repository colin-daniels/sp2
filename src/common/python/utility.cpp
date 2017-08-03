#include <Python.h>

#include "utility.hpp"
#include "error.hpp"
#include "py_scoped_t.hpp"
#include "py_opaque_t.hpp"

namespace sp2 {
namespace python {

py_opaque_t merge_dictionaries(const py_opaque_t &a,
    const py_opaque_t &b, merge_strategy strategy)
{
    auto size = [](auto &dict) {
        Py_ssize_t s = PyDict_Size(dict.inner().raw());
        throw_on_py_err("error obtaining Python dictionary size");
        return s;
    };

    auto merge_into_copy = [&](bool override) {
        auto c = scope(PyDict_Copy(a.inner().raw()));
        throw_on_py_err("error copying Python dictionary");

        PyDict_Merge(c.raw(), b.inner().raw(), override ? 1:0);
        throw_on_py_err("error merging Python dictionaries");

        return c;
    };

    switch (strategy) {
    case merge_strategy::USE_FIRST:  return opaque(merge_into_copy(false));
    case merge_strategy::USE_SECOND: return opaque(merge_into_copy(true));
    case merge_strategy::ERROR:
        auto c = merge_dictionaries(a, b, merge_strategy::USE_FIRST);

        if (size(c) != size(a) + size(b))
            throw std::runtime_error("Conflicting key in dict merge!");

        return c;
    }
}

} // namespace python
} // namespace sp2