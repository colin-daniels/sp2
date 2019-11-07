#include <Python.h> // Must be first include

#include "py_ref_t.hpp"
#include "py_object_t.hpp"
#include "common/python/error.hpp"

using namespace std;

namespace sp2 {
namespace python {

py_ref_t::py_ref_t(PyObject *o)
    : obj(o)
{}

py_ref_t::operator bool() const
{
    return (bool) obj;
}

py_ref_t::py_ref_t(py_ref_t &&other)
    : obj(other.steal())
{}

py_ref_t &py_ref_t::operator=(py_ref_t &&other)
{
    if (this != &other)
    {
        if (obj)
        {
            // we could implicitly destroy the existing reference, but with no
            // clear use case, the conservative choice is to require explicit
            // destruction
            throw logic_error("attempted to overwrite occupied py_ref_t");
        }
        obj = other.steal();
    }
    return *this;
}

py_ref_t::~py_ref_t()
{
    destroy();
}

py_ref_t py_ref_t::dup() const
{
    Py_XINCREF(obj);
    return py_ref_t(obj);
}

PyObject *py_ref_t::raw() const
{
    return obj;
}

py_ref_t py_ref_t::move()
{
    return std::move(*this);
}

PyObject *py_ref_t::steal()
{
    auto tmp = obj;
    obj = NULL;
    return tmp;
}

void py_ref_t::destroy()
{
    Py_XDECREF(obj);
    obj = NULL;
}

py_ref_t &py_ref_t::operator=(const py_ref_t &other)
{
    *this = other.dup();
    return *this;
}

py_ref_t::py_ref_t(const py_ref_t &other)
    :py_ref_t(other.dup())
{}

py_ref_t scope(PyObject *o)
{
    return py_ref_t(o);
}

py_ref_t scope_dup(PyObject *o)
{
    Py_XINCREF(o);
    return py_ref_t(o);
}

// --------------------------------

// Implementation of repr() and str().
// F is a function(PyObject *) -> PyObject * returning a new reference to a unicode 'str' object
template<typename F>
std::string str_impl(const py_ref_t &o, F stringify)
{

    auto py_str = scope(stringify(o.raw()));
    throw_on_py_err("repr: error stringifying");

    // NOTE: we are not responsible for deallocating this;
    //       it's lifetime is bound to the py_str python object.
    const char *str = PyUnicode_AsUTF8(py_str.raw());
    throw_on_py_err("repr: error encoding as utf8");

    return string(str);
}

string str(const py_ref_t &o)
{
    return str_impl(o, [](auto x) { return PyObject_Str(x); });
}

string repr(const py_ref_t &o)
{
    return str_impl(o, [](auto x) { return PyObject_Repr(x); });
}

py_ref_t getattr(const py_ref_t &o, const char *attr)
{
    auto tmp = scope(PyObject_GetAttrString(o.raw(), attr));
    throw_on_py_err();
    return tmp;
}

py_ref_t getattr(const py_ref_t &o, const char *attr, const py_ref_t &def)
{
    if (hasattr(o, attr))
        return getattr(o, attr);
    else
        return def;
}

bool hasattr(const py_ref_t &o, const char *attr)
{
    // From the python docs, "this always succeeds"
    return bool(PyObject_HasAttrString(o.raw(), attr));
}

void setattr(py_ref_t &o, const char *attr, const py_ref_t &value)
{
    PyObject_SetAttrString(o.raw(), attr, value.raw());
    throw_on_py_err(("error setting attribute '"s + attr + "'").c_str());
}

py_object_t opaque(py_ref_t &&scoped)
{
    typedef py_object_t::impl_t impl_t;
    return py_object_t{std::make_shared<impl_t>(std::move(scoped))};
}

py_object_t opaque(py_ref_t &scoped)
{
    return opaque(scoped.dup());
}

} // namespace python
} // namespace sp2
