#include <Python.h> // Must be first include

#include "py_opaque_t.hpp"
#include "py_scoped_t.hpp"
#include "common/python/error.hpp"

#include "common/python/modules/fake_modules.hpp"

using namespace sp2;
using sp2::python::py_opaque_t;
using sp2::python::py_scoped_t;

// anonymous namespace for private things
namespace {

/* --------------------------------------------------------------------- */
// helper type for python function arguments

// a pair of *args and **kw
typedef std::pair<py_scoped_t, py_scoped_t> args_t;

args_t args_and_kw(py_scoped_t args, py_scoped_t kw)
{
    return std::make_pair(std::move(args), std::move(kw));
}

args_t just_args(py_scoped_t args)
{
    auto kw = python::scope(PyDict_New());
    return std::make_pair(std::move(args), std::move(kw));
}

args_t just_kw(py_scoped_t kw)
{
    auto args = python::scope(PyTuple_New(0));
    return std::make_pair(std::move(args), std::move(kw));
}

/* --------------------------------------------------------------------- */
// helper functions

/// Invoke an object's __call__.
py_scoped_t call_callable(py_scoped_t &function, const args_t &argpair)
{
    auto &args = argpair.first;
    auto &kw = argpair.second;
    if (!(args && kw))
        throw std::logic_error(
            "tried to call function with invalid or previously used args_t");

    if (!PyArg_ValidateKeywordArguments(kw.raw()))
        throw std::logic_error("invalid keyword arguments");
    python::throw_on_py_err(); // errors not covered by the above

    if (!PyCallable_Check(function.raw()))
        throw std::logic_error("not callable");

    auto retval = python::scope(
        PyObject_Call(function.raw(), args.raw(), kw.raw()));

    python::throw_on_py_err("error calling python callable");
    return retval;
}

// Overload for convenience when there are both kw and args.
// (having only one would be ambiguous)
py_scoped_t call_callable(py_scoped_t &function, const py_scoped_t &args,
    const py_scoped_t &kw)
{
    auto argpair = args_and_kw(args.dup(), kw.dup());
    return call_callable(function, move(argpair));
}

py_scoped_t call_callable_nullable(py_scoped_t &function,
    const py_scoped_t &args,
    const py_scoped_t &kw)
{
    auto true_args = args;
    auto true_kw = kw;
    if (!true_args)
        true_args = python::scope(Py_BuildValue("()"));
    if (!true_kw)
        true_kw = python::scope(Py_BuildValue("{}"));
    return call_callable(function, args_and_kw(true_args, true_kw));
}

} // anonymous namespace

py_scoped_t& sp2::python::py_opaque_t::inner()
{
    if (!_impl)
    {
        throw std::logic_error("detected py_opaque_t with a NULL unique_ptr!"
            " The inner py_scoped_t should be NULL instead.");
    }

    return *_impl;
}

const py_scoped_t& sp2::python::py_opaque_t::inner() const
{
    if (!_impl)
    {
        throw std::logic_error("detected py_opaque_t with a NULL unique_ptr!"
            " The inner py_scoped_t should be NULL instead.");
    }

    return *_impl;
}

py_opaque_t::py_opaque_t()
    // make a wrapped nullptr, not an actual nullptr
    : _impl(std::make_shared<py_scoped_t>(py_scoped_t()))
{ }

py_opaque_t::~py_opaque_t() = default;
py_opaque_t::py_opaque_t(py_opaque_t &&) noexcept = default;
py_opaque_t& py_opaque_t::operator=(py_opaque_t &&) noexcept = default;
py_opaque_t::py_opaque_t(const py_opaque_t &) = default;
py_opaque_t &py_opaque_t::operator=(const py_opaque_t &) = default;
py_opaque_t::py_opaque_t(std::shared_ptr<py_opaque_t::impl_t> &&impl)
    : _impl(std::move(impl))
{
    if (!_impl)
    {
        throw std::logic_error("cannot construct py_opaque_t from "
            "NULL unique_ptr. (the py_scoped_t should be null instead)");
    }
}

bool py_opaque_t::hasattr(const char *attr) const
{ return python::hasattr(inner(), attr); }

py_opaque_t py_opaque_t::getattr(const char *attr) const
{ return opaque(python::getattr(inner(), attr)); }

py_opaque_t py_opaque_t::getattr(const char *attr, const py_opaque_t &def) const
{ return opaque(python::getattr(inner(), attr, def.inner())); }

void py_opaque_t::setattr(const char *attr, const py_opaque_t &value)
{ return python::setattr(inner(), attr, value.inner()); }

py_opaque_t py_opaque_t::call(const py_opaque_t &args, const py_opaque_t &kw)
{ return opaque(call_callable_nullable(inner(), args.inner(), kw.inner())); }

void py_opaque_t::destroy()
{ return inner().destroy(); }

py_opaque_t::operator bool() const
{ return bool(inner()); }

py_opaque_t py_opaque_t::_tuple(const std::vector<py_opaque_t> &ts)
{
    // IMPORTANT:
    //
    // This implementation MUST NOT use py_opaque_t conversion methods,
    //  because it is used extensively in their implementation!

    std::vector<py_scoped_t> arg_inners;
    for (const py_opaque_t &arg : ts)
        arg_inners.push_back(arg.inner());

    // Make a list.
    //
    // This invokes the specific overload:
    //     to_python(const vector<py_scoped_t>&, py_scoped_t&)
    // which can be manually verified to never involve a py_opaque_t conversion.
    py_scoped_t py_args_list;
    if (!to_python(arg_inners, py_args_list))
        throw std::runtime_error("Error creating argument list.");

    // Make a tuple
    auto py_args = scope(PySequence_Tuple(py_args_list.raw()));
    throw_on_py_err();
    return opaque(py_args);
}

py_opaque_t sp2::python::import(const char *mod_name)
{
    auto name = scope(PyUnicode_DecodeFSDefault(mod_name));
    throw_on_py_err();

    auto module = scope(PyImport_Import(name.raw()));
    throw_on_py_err();
    return opaque(module);
}

py_opaque_t sp2::python::import(const std::string &mod_name)
{ return import(mod_name.c_str()); }


