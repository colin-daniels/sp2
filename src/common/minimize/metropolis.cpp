#include "common/python.hpp"

namespace sp2 {
namespace minimize {
namespace metropolis {

#ifdef SP2_ENABLE_PYTHON

using namespace std;
using namespace sp2::python;

sp2::python::py_object_t make_param_pack(vector<double> carts,
    const double lattice[3][3], vector<double> force)
{

    // interpret as 3N cartesian coords
    size_t width = 3;
    if (carts.size() % width != 0)
        throw logic_error("vector not divisible by width");

    size_t height = carts.size() / width;

    auto py_carts = py_from(as_ndarray(carts, {height, width}));

    auto c_lattice = vector<double>(&lattice[0][0], &lattice[0][0] + 9);
    auto py_lattice = py_from(as_ndarray(c_lattice, {3, 3}));

    auto py_force = py_from(as_ndarray(force, {height, width}));

    auto kw = py_import("builtins").getattr("dict").call();
    auto setter = kw.getattr("__setitem__");
    setter.call(py_tuple(py_from("carts"s), py_carts), {});
    setter.call(py_tuple(py_from("lattice"s), py_lattice), {});
    setter.call(py_tuple(py_from("force"s), py_force), {});

    return kw;
}

#endif // SP2_ENABLE_PYTHON

} // namespace metropolis
} // namespace minimize
} // namespace sp2