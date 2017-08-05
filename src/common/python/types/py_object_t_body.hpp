#ifndef SP2_PY_OBJECT_T_BODY_HPP
#define SP2_PY_OBJECT_T_BODY_HPP

#include "py_object_t_fwd.hpp"

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

// If an incomplete type will suffice, consider including the FWD instead.

#include "py_ref_t_fwd.hpp"

#include <vector>
#include <string>
#include <memory>

namespace sp2 {
namespace python {

/// High-level wrapper around python objects.
///
/// This type can be used to perform operations on python objects without
/// needing to import Python headers.
class py_object_t
{
public:
    /// pImpl target type.
    typedef py_ref_t impl_t;

private:
    std::shared_ptr<impl_t> _impl;

public:
    py_ref_t& inner();
    const py_ref_t& inner() const;

    py_object_t(std::shared_ptr<py_ref_t>&&);

    py_object_t();
    py_object_t(py_object_t&&) noexcept;
    py_object_t& operator=(py_object_t&&) noexcept;
    ~py_object_t();
    py_object_t(const py_object_t&);
    py_object_t& operator=(const py_object_t&);

    explicit operator bool() const;

    void destroy();

    //-----------------------
    // Provide a core of high-level functionality.
    //
    // Together, attribute access and the ability to call functions enables 99%
    // of all use cases for manipulating python objects and lets us move a lot
    // of logic out of bindings.cpp and closer to where it belongs.

    /// Test if a python object has an attribute.
    ///
    /// Equivalent to 'hasattr(obj, name)'. Never fails.
    // FIXME though I don't know what it does if your object is NULL...
    bool hasattr(const char *attr) const;
    bool hasattr(const std::string &attr) const
    { return hasattr(attr.c_str()); }

    /// Access an attribute of a python object.
    ///
    /// Equivalent to 'getattr(obj, name)'.
    /// Throws an exception if the attribute does not exist.
    py_object_t getattr(const char *attr) const;
    py_object_t getattr(const std::string &attr) const
    { return getattr(attr.c_str()); }

    /// Access an attribute of a python object, or a default value.
    ///
    /// Equivalent to 'getattr(obj, name, def)'.
    py_object_t getattr(const char *attr, const py_object_t &def) const;
    py_object_t getattr(const std::string &attr, const py_object_t &def) const
    { return getattr(attr.c_str(), def); }

    /// Set an attribute of a python object.
    void setattr(const char *attr, const py_object_t &value);
    void setattr(const std::string &attr, const py_object_t &value)
    { return setattr(attr.c_str(), value); }

    /// Call a python callable.
    ///
    /// This takes an args tuple and a keyword dict. Either or both can be null,
    /// unlike *cough* some APIs...
    py_object_t call(const py_object_t &args, const py_object_t &kw);

    py_object_t call()
    { return call(py_object_t{}, {}); }

    /// Call a python callable.
    // This overload is not just for convenience; generic conversions to/from
    // py_object_t have no recourse but to call python functions, so this
    // overload solves the chicken-and-egg problem of building an args tuple.
    py_object_t call(std::vector<py_object_t> args)
    { return call(detail::py_tuple_noconv(args), {}); }

    /// Value-returning wrapper around from_python.
    /// This will require an explicit type annotation.
    ///
    /// Errors are communicated by std::runtime_error.
    template<
        typename T,
        typename = std::enable_if_t<std::is_default_constructible<T>::value>
    >
    T parse_as(const char *msg) const;

    template<
        typename T,
        typename = std::enable_if_t<std::is_default_constructible<T>::value>
    >
    T parse_as() const
    { return parse_as<T>("an error occurred converting data from python"); }

private:
    void debug_null() const;
};

} // namespace python
} // namespace sp2

#endif // header guard