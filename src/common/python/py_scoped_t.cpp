#include <Python.h> // Must be first include

#include "py_scoped_t.hpp"
#include "py_opaque_t.hpp"
#include "error.hpp"

using namespace std;

namespace sp2 {
namespace python {

py_scoped_t::py_scoped_t(PyObject *o)
    : obj(o)
{}

py_scoped_t::operator bool() const
{
    return (bool) obj;
}

py_scoped_t::py_scoped_t(py_scoped_t &&other)
    : obj(other.steal())
{}

py_scoped_t &py_scoped_t::operator=(py_scoped_t &&other)
{
    if (this != &other)
    {
        if (obj)
        {
            // we could implicitly destroy the existing reference, but with no
            // clear use case, the conservative choice is to require explicit
            // destruction
            throw logic_error("attempted to overwrite occupied py_scoped_t");
        }
        obj = other.steal();
    }
    return *this;
}

py_scoped_t::~py_scoped_t()
{
    destroy();
}

py_scoped_t py_scoped_t::dup() const
{
    Py_XINCREF(obj);
    return py_scoped_t(obj);
}

PyObject *py_scoped_t::raw() const
{
    return obj;
}

py_scoped_t py_scoped_t::move()
{
    return std::move(*this);
}

PyObject *py_scoped_t::steal()
{
    auto tmp = obj;
    obj = NULL;
    return tmp;
}

void py_scoped_t::destroy()
{
    Py_XDECREF(obj);
    obj = NULL;
}

py_scoped_t &py_scoped_t::operator=(const py_scoped_t &other)
{
    *this = other.dup();
}

py_scoped_t::py_scoped_t(const py_scoped_t &other)
    :py_scoped_t(other.dup())
{}

py_scoped_t scope(PyObject *o)
{
    return py_scoped_t(o);
}

py_scoped_t scope_dup(PyObject *o)
{
    Py_XINCREF(o);
    return py_scoped_t(o);
}

// --------------------------------

// Implementation of repr() and str().
// F is a function(PyObject *) -> PyObject * returning a new reference to a unicode 'str' object
template<typename F>
std::string str_impl(const py_scoped_t &o, F stringify)
{

    auto py_str = scope(stringify(o.raw()));
    throw_on_py_err("repr: error stringifying");

    // NOTE: we are not responsible for deallocating this;
    //       it's lifetime is bound to the py_str python object.
    char *str = PyUnicode_AsUTF8(py_str.raw());
    throw_on_py_err("repr: error encoding as utf8");

    return string(str);
}

string str(const py_scoped_t &o)
{
    return str_impl(o, [](auto x) { return PyObject_Str(x); });
}

string repr(const py_scoped_t &o)
{
    return str_impl(o, [](auto x) { return PyObject_Repr(x); });
}

py_scoped_t getattr(const py_scoped_t &o, const char *attr)
{
    auto tmp = scope(PyObject_GetAttrString(o.raw(), attr));
    throw_on_py_err();
    return move(tmp);
}

py_scoped_t getattr(const py_scoped_t &o, const char *attr, const py_scoped_t &def)
{
    if (hasattr(o, attr))
        return getattr(o, attr);
    else
        return move(def);
}

bool hasattr(const py_scoped_t &o, const char *attr)
{
    // From the python docs, "this always succeeds"
    return bool(PyObject_HasAttrString(o.raw(), attr));
}

void setattr(py_scoped_t &o, const char *attr, const py_scoped_t &value)
{
    PyObject_SetAttrString(o.raw(), attr, value.raw());
    throw_on_py_err(("error setting attribute '"s + attr + "'").c_str());
}

py_opaque_t opaque(py_scoped_t &&scoped)
{
    typedef py_opaque_t::impl_t impl_t;
    return py_opaque_t{std::make_shared<impl_t>(std::move(scoped))};
}

py_opaque_t opaque(py_scoped_t &scoped)
{
    return opaque(scoped.dup());
}

} // namespace python
} // namespace sp2
