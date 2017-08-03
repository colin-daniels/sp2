#include "Python.h" // must be first include

#include "base_generic_raw.hpp"
#include "common/python/py_scoped_t.spp"
#include "common/python/bindings.hpp" // make_param_pack

using namespace std;

namespace sp2 {
namespace python {

/* -------------------------------------------------------------------------- */
// Helper implementation.
// Serializes to/from string.
template<typename E>
bool from_python_by_enum_map(const py_scoped_t &py, E &out, E null_value,
    const char *message)
{
    std::string s;
    if (!from_python(py, s))
        return false;

    out = enum_from_str(s, null_value);
    if (out == null_value)
    {
        std::cerr << message << std::endl;
        return false;
    }
    return true;
}

/* -------------------------------------------------------------------------- */

bool from_python_concrete(const py_scoped_t &py, structural_mutation_type &c)
{
    return from_python_by_enum_map(py, c,
        structural_mutation_type::INVALID,
        (std::string()
         + "Invalid value for structural mutation type: "
         + repr(py)
        ).c_str());
}

bool from_python_concrete(const py_scoped_t &py, structural_mutation_t &c)
{
    tuple<structural_mutation_type, as_ndarray_t<double>> tup;
    if (!from_python(py, tup))
        return false;

    auto type = get<0>(tup);
    auto data = get<1>(tup);

    c = {type, data};
    return true;
}

} // namespace python
} // namespace sp2
