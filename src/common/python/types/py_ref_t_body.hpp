#ifndef SP2_PY_REF_T_BODY_HPP
#define SP2_PY_REF_T_BODY_HPP

#include "py_ref_t_fwd.hpp"

// NOTE: This is a TYPE-BODY-ONLY header file.
//
// It only defines the following:
//  * The type.
//  * Member functions which can be defined entirely in terms of other
//    member functions. (plus light dependencies that have no risk of cycles)
//
// Consumers within the library may prefer to include this if they require a
// complete type but do not need to use template member functions or any of the
// free functions provided in the HPP file. This reduces the #include footprint,
// which can help eliminate some dependency cycles.

// If an incomplete type will suffice, consider including the FWD instead;
// it has the advantage of not requiring Python headers.

#include "diagnostic/expect_python_headers"

namespace sp2 {
namespace python {

/// The designated owner of a Python object reference.
///
/// This type is for internal use in parts of the codebase that interact
/// directly with the CPython API.  It is a scoped reference to a python object
/// that uses RAII to handle Py_DECREF.
///
/// This makes it somewhat easier to reason about exception safety,
/// though it is not a panacea.
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
/// 'const py_ref' won't allow you to modify its pointer to point
/// somewhere else, but you are free to modify its referent in any other way.
/// This admission is made because a 'const PyObject*' is nearly useless.
///
/// The contained object may be NULL.
class py_ref_t
{
    mutable PyObject *obj = nullptr;

public:

    /// null constructor
    py_ref_t() {};

    /// PyObject constructor
    ///
    /// WARNING: The current implementation assumes it is being given ownership
    ///           of the pointer, because this is the least expensive default
    ///           and nobody is expected to be using this.
    ///          Unfortunately this is hostile to modern C++ convention.
    ///
    /// FIXME: Actually, I can't even recall why this is public.
    ///        I know there was a reason, but IIRC it was just an
    ///         implementation detail, so it should be documented...
    explicit py_ref_t(PyObject *o);

    explicit operator bool() const;

    py_ref_t &operator=(const py_ref_t &other);
    py_ref_t(const py_ref_t &other);

    /// move constructor
    py_ref_t(py_ref_t &&other);

    /// move assignment operator
    py_ref_t& operator=(py_ref_t &&other);

    ~py_ref_t();

    /// Increment the refcount and return a new scoped reference.
    /// FIXME this is pointless now that there's a copy constructor
    py_ref_t dup() const;

    /// Borrow the reference without touching the refcount.
    ///
    /// This is the appropriate method for interfacing with most Python APIs.
    ///
    /// For reasons discussed earlier, the returned pointer is mutable.
    /// Just... be reasonable, okay?
    PyObject *raw() const;

    /// Leak the reference, preventing the DECREF that would otherwise occur
    /// at scope exit. The scoped reference will become NULL.
    ///
    /// Necessary for working with Python API functions that steal references,
    /// such as PyTuple_SetItem.
    PyObject *steal();

    /// Postfix move() for convenience, and for easier interchange with dup().
    ///
    /// Note a small semantic difference in that, due to the by-value return
    /// type, this is guaranteed to clear the receiver.
    py_ref_t move();

    /// Explicit destructor.
    ///
    /// Destroy the reference early, decrementing the refcount and nulling out
    /// the pointer so that nothing happens at scope exit.
    ///
    /// This can be used to explicitly control the destruction order in
    /// places where the natural order of destruction would not be safe.
    void destroy();
};

} // namespace python
} // namespace sp2
#endif // header guard