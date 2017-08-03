/* ========================================================================== */
#ifndef SP2_PY_SCOPED_T_TPP
#define SP2_PY_SCOPED_T_TPP

// no forward dependencies

#undef SP2_PY_SCOPED_T_TPP
#else
#include "diagnostic/forward_dependency_cycle"
#endif
/* ========================================================================== */

namespace sp2 {
namespace python {

/// The designated owner of a Python object reference.
///
/// This scoped reference to a python object that uses RAII to handle Py_DECREF.
/// This makes it somewhat easier to reason about exception safety, though it is
/// not a panacea.
///
/// One should still be careful to consider destruction order (since a decref
/// can potentially invoke arbitrary python code), and read the Python API docs
/// carefully to understand when references are duplicated, borrowed, and
/// stolen.
///
/// The default PyObject * constructor does NOT perform an incref, since the
/// majority of Python API functions return a freshly-incremented reference.
/// For those rare functions that return borrowed references, you should
/// use the explicit 'scope_dup' constructor instead.
///
/// 'const' guarantees for this type are particularly weak. In general, a
/// 'const py_scoped_t' won't allow you to modify its pointer to point
/// somewhere else, but you are free to modify its referent in any other way.
/// This admission is made because a 'const PyObject*' is nearly useless.
///
/// The contained object may be NULL.
class py_scoped_t;

} // namespace python
} // namespace sp2
