/* ========================================================================== */
#ifndef SP2_PY_OBJECT_T_FWD_DAG
#define SP2_PY_OBJECT_T_FWD_DAG

// no forward dependencies

#undef SP2_PY_OBJECT_T_FWD_DAG
#else
#include "diagnostic/forward_dependency_cycle"
#endif
/* ========================================================================== */

#include <vector>
#include <string>

namespace sp2 {
namespace python {

class py_object_t;

/// Import a module as a python object.
///
/// This is one of the fundamental methods of obtaining a 'py_object'.
/// (the other being to convert data from C++ via 'py_from' or 'py_tuple').
///
/// This imports a module by its (package-qualified) name through the standard
/// python module-loading machinery, which caches imports. In other words,
/// this MAY invoke module.__init__() (and therefore one should be prepared for
/// arbitrary python errors), but it will only do so once.
py_object_t py_import(const char *mod_name);
py_object_t py_import(const std::string &mod_name);

namespace detail {

// a monomorphic variant of py_tuple that is guaranteed not to perform
// conversions between py_object_t and arbitrary types;
py_object_t py_tuple_noconv(const std::vector<py_object_t> &ts);

} // namespace detail
} // namespace python
} // namespace sp2