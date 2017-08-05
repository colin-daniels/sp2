#ifndef SP2_PY_REF_T_HPP
#define SP2_PY_REF_T_HPP

#include "py_ref_t_body.hpp"

#include "py_object_t_fwd.hpp"

#include "diagnostic/expect_python_headers"

#include <string>
#include <utility>
#include <stdexcept>

namespace sp2 {
namespace python {

/// explicit constructor from a new ref
py_ref_t scope(PyObject *o);

/// explicit constructor from a borrowed ref, which makes a new reference
py_ref_t scope_dup(PyObject *o);

// --------------------------------
// less typing to construct 'py_object_t's

py_object_t opaque(py_ref_t &&scoped);

py_object_t opaque(py_ref_t &scoped);

// --------------------------------

/// get an object's repr() in utf8, mostly for debug purposes.
std::string repr(const py_ref_t &o);

/// get an object's str() in utf8, mostly for debug purposes.
std::string str(const py_ref_t &o);

// --------------------------------

/// Test if a python object has an attribute.
///
/// Equivalent to 'hasattr(obj, name)'.  Never fails.
bool hasattr(const py_ref_t &o, const char *attr);

/// Access an attribute of a python object.
///
/// Equivalent to 'getattr(obj, name)'.
/// Throws an exception if the attribute does not exist.
py_ref_t getattr(const py_ref_t &o, const char *attr);

/// Access an attribute of a python object, or a default value.
///
/// Equivalent to 'getattr(obj, name, def)'.
/// Throws an exception if the attribute does not exist.
py_ref_t getattr(const py_ref_t &o, const char *attr, const py_ref_t &def);

/// Set an attribute of a python object.
void setattr(py_ref_t &o, const char *attr, const py_ref_t &def);

} // namespace python
} // namespace sp2

#endif // header guard
