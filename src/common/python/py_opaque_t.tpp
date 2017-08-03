/* ========================================================================== */
#ifndef SP2_PYTHON_PY_OPAQUE_T_TPP
#define SP2_PYTHON_PY_OPAQUE_T_TPP

// no forward dependencies

#undef SP2_PYTHON_PY_OPAQUE_T_TPP
#else
#include "diagnostic/forward_dependency_cycle"
#endif
/* ========================================================================== */

#include <vector>
#include <string>

namespace sp2 {
namespace python {

/// An opaque, pImpl-style wrapper around a python pointer
/// for safe encapsulation even in code that does not import
/// the cPython headers.
class py_opaque_t;

/// Import a module as a python object.
///
/// This is one of the fundamental methods of obtaining a 'py_opaque_t'.
/// (the other being to convert values from C++ type).
///
/// This imports a module by its (package-qualified) name through the standard
/// python module-loading machinery, which caches imports. In other words,
/// this MAY invoke module.__init__() (and therefore one should be prepared for
/// arbitrary python errors), but it will only do so once.
py_opaque_t import(const char *mod_name);
py_opaque_t import(const std::string &mod_name);

} // namespace python
} // namespace sp2