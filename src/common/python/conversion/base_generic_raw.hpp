#ifndef SP2_PYTHON_CONVERSION_GENERIC_RAW_HPP
#define SP2_PYTHON_CONVERSION_GENERIC_RAW_HPP

#include "base_generic_raw_fwd.hpp"

// base_generic_raw.hpp:
//
//  Provides '{to,from}_python' on generic types for 'py_ref_t'.
//  These implementations are able to directly use CPython API functions on the
//   corresponding containers in python.

//      !!!!!!!!!!!!!!!!!!!
//        !!  CAUTION  !!
//
// 'base_generic.hpp' and 'base_generic_raw.hpp' both independently implement
// the same functions on different types.  Never edit one file without editing
// the other.
//
// Unfortunately there is currently no viable abstraction that would enable
// a single implementation to be written for both types.
//
//        !!  CAUTION  !!
//      !!!!!!!!!!!!!!!!!!!

#include "diagnostic/expect_python_headers"

#include "common/python/error.hpp"
#include "common/python/types/py_ref_t.hpp"
#include "base_monomorphic.hpp"

#include <iostream>

namespace sp2 {
namespace python {

/* --------------------------------------------------------------------- */
// Generic conversions of:
//
//  - std::vector ---> list
//  - std::vector <--- sequence or iterable

template<typename T>
bool to_python(const std::vector<T> &vec, py_ref_t &list);

template<typename T>
bool to_python(const std::vector<T> &vec, py_ref_t &list)
{
    if (list)
        throw std::logic_error(
            "attempted to serialize into occupied py_ref_t");

    // of course, pre-allocating a list of the right size would be faster,
    // but 'new' and 'append' keeps the object in a consistent state
    list = scope(PyList_New(0));
    throw_on_py_err("unexpected error constructing python list");

    for (auto &x : vec)
    {
        py_ref_t px;
        if (!to_python(x, px))
        {
            // ...and here's why it's nice to keep 'list' in a consistent state
            list.destroy();
            return false;
        }

        PyList_Append(list.raw(), px.raw());
        if (print_on_py_err())
        {
            list.destroy();
            throw std::runtime_error(
                "unexpected error appending to python list");
        }
    }

    return true;
}

template<typename T>
bool from_python(const py_ref_t &o, std::vector<T> &vec)
{
    vec.clear();

    // NOTICE
    // Although the docs for PySequence_Fast say that it takes a sequence,
    // it actually accepts iterables as well.
    //
    // It's basically:  x if isinstance(x, (list, tuple)) else list(x)
    //
    // Both of these cases are fine; except that iterables include iterators,
    // which we'd rather NOT accept since they would get consumed the first time
    // we fail to deserialize them.
    if (PyIter_Check(o.raw()))
    {
        std::cerr << "error: refusing to deserialize a one-time use iterator. "
                  << "Please use e.g. list() to strictly evaluate it."
                  << std::endl;
        return false;
    }

    // With that out of the way, any errors from PySequence_Fast are
    // legitimate concerns.
    py_ref_t seq = scope(PySequence_Fast(o.raw(), "expected a sequence!"));
    throw_on_py_err();

    Py_ssize_t size = PySequence_Fast_GET_SIZE(seq.raw());
    vec.reserve(size);
    for (Py_ssize_t i = 0; i < size; i++)
    {
        // PySequence_Fast_GET_ITEM returns a borrowed reference,
        //  and does not ever (?) throw python exceptions.
        auto py_item = scope_dup(PySequence_Fast_GET_ITEM(seq.raw(), i));

        T item;
        if (!from_python(py_item, item))
            return false;

        vec.push_back(std::move(item));
    }

    return true;
}

/* --------------------------------------------------------------------- */
// Generic conversions of:
//
//  - std::tuple ---> list (FIXME: justify why not tuple)
//  - std::tuple <--- sequence or iterable

template<std::size_t Idx, class ...Ts>
constexpr bool _from_python_rec_raw(Ts &&...)
{ return true; }

template<std::size_t Idx = 0, class ...Ts,
    class = std::enable_if_t<(Idx < sizeof...(Ts))>>
bool _from_python_rec_raw(const std::vector<py_ref_t> &pys, std::tuple<Ts...> &cs)
{
    return from_python(pys[Idx], std::get<Idx>(cs)) &&
           _from_python_rec_raw<Idx + 1>(pys, cs);
}

template<typename... Ts>
bool from_python(const py_ref_t &py, std::tuple<Ts...> &cs)
{
    std::vector<py_ref_t> pys;
    if (!from_python(py, pys))
        return false;

    if (pys.size() != sizeof...(Ts))
    {
        std::cerr << "expected python sequence of length " << sizeof...(Ts)
                  << ", got length " << pys.size() << std::endl;
        return false;
    }

    return _from_python_rec_raw(pys, cs);
}

// -----

template<std::size_t Idx>
constexpr bool _to_python_rec_raw(...)
{ return true; }

template<std::size_t Idx = 0, class ...Ts,
    class = std::enable_if_t<(Idx < sizeof...(Ts))>>
bool _to_python_rec_raw(const std::tuple<Ts...> &cs, std::vector<py_ref_t> &pys)
{
    return to_python(std::get<Idx>(cs), pys[Idx]) &&
           _to_python_rec_raw<Idx + 1>(cs, pys);
}

template<typename... Ts>
bool to_python(const std::tuple<Ts...> &cs, py_ref_t &py)
{
    std::vector<py_ref_t> pys(sizeof...(Ts));

    if (!_to_python_rec_raw(cs, pys))
    {
        // destruct any python references that were created prior to the error
        pys.resize(0);
    }

    return to_python(pys, py);
}

/* --------------------------------------------------------------------- */
} // namespace python
} // namespace sp2

#endif // SP2_PYTHON_CONVERSION_GENERIC_RAW_HPP
