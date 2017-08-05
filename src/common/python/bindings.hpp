#ifndef SP2_PYTHON_BINDINGS_HPP
#define SP2_PYTHON_BINDINGS_HPP

#include "common/python/types/py_object_t_fwd.hpp"

// This header is an odd grabbag of logic for things outside the python
// library code, but whose implementations call into the CPython API.

namespace sp2 {
namespace python {

// temp
namespace structural_metropolis {

py_object_t make_param_pack(std::vector<double> carts,
    const double lattice[3][3], std::vector<double> force);

}

namespace run_phonopy {

#warning misplaced
/// Produce extra structural_metropolis kw args for run_phonopy.
///
/// These will not conflict with any of structural_metropolis' own kw args.
py_object_t make_extra_kw(std::vector<size_t> sc_to_prim);

} // namespace run_phonopy

} // namespace python
} // namespace sp2

#endif // SP2_PYTHON_BINDINGS_HPP
