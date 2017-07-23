#ifndef SP2_PYTHON_CONVERSION_HPP
#define SP2_PYTHON_CONVERSION_HPP

#ifndef Py_PYTHON_H
#error This module is designed for internal use by source files in common/python, \
and is not safe to use unless Python.h is already included before all standard \
library headers.
#endif // Py_PYTHON_H

#include "common/python/util.hpp"
#include "common/python/numpy_util.hpp"
#include "common/python/include_numpy.hpp"
#include "phonopy/structural_mutation.hpp"

#include <cstdint>
#include <vector>
#include <iostream>

// to_python and from_python.
//
// 'to_python' turns a variety of C++ data types into python equivalents.
// 'from_python' reads python data back into strongly typed C++.
//

// to_python:
//
// to_python generally tries to produce types which approximate or represents
// the C++ type in spirit.  For cases that are unclear, selection of the
// desired python type may be accomplished through the use of wrapper types.
// (though these don't exist yet)
//
// to_python is generally unlikely to fail, but when it does, all references
// created by it are guaranteed to be cleaned up, and it will leave the
// output argument NULL and return false.

// from_python:
//
// 'from_python' must be able to recover an object from the output of
// 'to_python', but beyond that it is fairly permissive in what it accepts,
// using Python APIs that provide the conversions typically accepted in
// python (e.g. accepting an int in place of a float, or an iterable in place
// of a list, etc.)
//
// On failure, from_python returns false, and leaves the output argument in
// an unspecified state.
//
// Failure is not uncommon, and in our case generally results from poor output
// returned by a user script; hence good error messages are a must.
// The current implementations typically try to take advantage of python's
//  own fantastic error messages and tracebacks by using PyErr_Print
//  before returning false.
//
// Of course, unconditionally printing errors makes the current implementations
// unsuitable for use cases such as "attempt to deserialize as type A, then
// attempt to deserialize as type B if that fails". This may be addressed in
// the future.

namespace sp2 {
namespace python {

bool to_python(const long &c, py_scoped_t &py);

bool to_python(const long long &c, py_scoped_t &py);

bool to_python(const unsigned long &c, py_scoped_t &py);

bool to_python(const unsigned long long &c, py_scoped_t &py);

bool to_python(const double &c, py_scoped_t &py);

bool to_python(const bool &c, py_scoped_t &py);

bool to_python(const std::string &c, py_scoped_t &py);

bool to_python(const std::nullptr_t &c, py_scoped_t &py);

//bool to_python(const structural_mutation_t &c, py_scoped_t &py);
bool to_python(py_scoped_t &c, py_scoped_t &py);

bool from_python(py_scoped_t &py, long &c);

bool from_python(py_scoped_t &py, long long &c);

bool from_python(py_scoped_t &py, unsigned long &c);

bool from_python(py_scoped_t &py, unsigned long long &c);

bool from_python(py_scoped_t &py, double &c);

bool from_python(py_scoped_t &py, bool &c);

bool from_python(py_scoped_t &py, std::string &c);

bool from_python(py_scoped_t &py, std::nullptr_t &c);

bool from_python(py_scoped_t &py, structural_mutation_t &c);

bool from_python(py_scoped_t &py, py_scoped_t &c);

/* --------------------------------------------------------------------- */
// Value returning conversions: {to,from}_python_strict
//
// These can be more convenient due to not taking an in/out parameter.
// They wrap {to,from}_python, turning conversion failures into exceptions.

/// Value-returning wrapper around to_python.
///
/// Errors are communicated by std::runtime_error.
template<typename T>
py_scoped_t to_python_strict(const T &c, const char *msg)
{
    py_scoped_t py;
    if (!to_python(c, py))
        throw std::runtime_error(msg);
    return py;
}

template<typename T>
py_scoped_t to_python_strict(const T &c)
{
    return to_python_strict(c,
        "an error occurred converting data into python objects");
}

/// Value-returning wrapper around from_python.
/// This will require explicit type annotations.
///
/// Errors are communicated by std::runtime_error.
template<
    typename T,
    typename = std::enable_if_t<std::is_default_constructible<T>::value>
>
T from_python_strict(py_scoped_t &py, const char *msg)
{
    auto c = T();
    if (!from_python(py, c))
        throw std::runtime_error(msg);
    return c;
}

template<
    typename T,
    typename = std::enable_if_t<std::is_default_constructible<T>::value>
>
T from_python_strict(py_scoped_t &py)
{
    return from_python_strict<T>(py,
        "an error occurred converting data from python");
}

/* --------------------------------------------------------------------- */
// Generic conversions of:
//
//  - std::vector ---> list
//  - std::vector <--- sequence or iterable

template<typename T>
bool to_python(const std::vector<T> &vec, py_scoped_t &list)
{
    if (list)
        throw std::logic_error(
            "attempted to serialize into occupied py_scoped_t");

    // of course, pre-allocating a list of the right size would be faster,
    // but 'new' and 'append' keeps the object in a consistent state
    list = scope(PyList_New(0));
    throw_on_py_err("unexpected error constructing python list");

    for (auto &x : vec)
    {
        py_scoped_t px;
        if (!to_python(x, px))
        {
            // ...and here's why it's nice to keep 'list' in a consistent state
            list.destroy();
            return false;
        }

        PyList_Append(list.raw(), px.raw());
        //std::cerr << "ITEMS: " << std::endl;
        //int n = PyList_GET_SIZE(list.raw());
        //for (int i=0; i < n; i++) {
        //    std::cerr << "- " << PyList_GET_ITEM(list.raw(), i) << std::endl;
        //}
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
bool from_python(py_scoped_t &o, std::vector<T> &vec)
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
    py_scoped_t seq = scope(PySequence_Fast(o.raw(), "expected a sequence!"));
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
// Generic conversions of as_ndarray_t <-> ndarray

/// Helper type to obtain the numpy dtype integral constant associated
/// with a data type.
template<typename T>
struct numpy_dtype {};

#define SPECIALIZE_NUMPY_DTYPE(T, VAL) \
template<> struct numpy_dtype<T> : std::integral_constant<int, VAL> {};

SPECIALIZE_NUMPY_DTYPE(bool, NPY_BOOL)
SPECIALIZE_NUMPY_DTYPE(int8_t, NPY_INT8)
SPECIALIZE_NUMPY_DTYPE(int16_t, NPY_INT16)
SPECIALIZE_NUMPY_DTYPE(int32_t, NPY_INT32)
SPECIALIZE_NUMPY_DTYPE(int64_t, NPY_INT64)
SPECIALIZE_NUMPY_DTYPE(uint8_t, NPY_UINT8)
SPECIALIZE_NUMPY_DTYPE(uint16_t, NPY_UINT16)
SPECIALIZE_NUMPY_DTYPE(uint32_t, NPY_UINT32)
SPECIALIZE_NUMPY_DTYPE(uint64_t, NPY_UINT64)
SPECIALIZE_NUMPY_DTYPE(float, NPY_FLOAT)
SPECIALIZE_NUMPY_DTYPE(double, NPY_DOUBLE)

// NOTE: DTYPE doubles as SFINAE
template<typename T, int DTYPE = numpy_dtype<T>::value>
bool to_python(const as_ndarray_t<T> &c, py_scoped_t &py)
{

    // copy data into a brand new array object.
    std::vector<npy_intp> shape(c.shape().begin(), c.shape().end());
    auto arr = scope(PyArray_SimpleNew(c.ndim(), shape.data(), DTYPE));
    if (print_on_py_err())
        return false;

    auto arr_data = (T *) PyArray_DATA((PyArrayObject *) arr.raw());
    throw_on_py_err("error accessing numpy array data");
    std::copy(c.data().begin(), c.data().end(), arr_data);

    py = std::move(arr);
    return true;
}

// NOTE: DTYPE doubles as SFINAE
template<typename T, int DTYPE = numpy_dtype<T>::value>
bool from_python(py_scoped_t &py, as_ndarray_t<T> &c)
{
    // Force the array into a contiguous layout if it isn't.
    int min_depth = 0; // ignore
    int max_depth = 0; // ignore
    auto contiguous = scope(
        PyArray_ContiguousFromAny(py.raw(), DTYPE, min_depth, max_depth));
    if (print_on_py_err())
        return false;

    auto arr = (PyArrayObject *) contiguous.raw();
    // NOTE: none of these set the python error state
    auto arr_data = (T *) PyArray_DATA(arr);
    auto arr_size = size_t(PyArray_SIZE(arr));
    auto arr_ndim = size_t(PyArray_NDIM(arr));
    auto *arr_dims = PyArray_DIMS(arr);

    std::vector<T> data(arr_data, arr_data + arr_size);
    std::vector<size_t> shape(arr_dims, arr_dims + arr_ndim);
    c = as_ndarray_t<T>(data, shape);
    return true;
}

/* --------------------------------------------------------------------- */
// Generic conversions of:
//
//  - std::tuple ---> list (FIXME: not tuple?)
//  - std::tuple <--- sequence or iterable

template<std::size_t Idx, class ...Ts>
constexpr bool _from_python_rec(Ts &&...)
{ return true; }

template<std::size_t Idx = 0, class ...Ts,
    class = std::enable_if_t<(Idx < sizeof...(Ts))>>
bool _from_python_rec(std::vector<py_scoped_t> &pys, std::tuple<Ts...> &cs)
{
    return from_python(pys[Idx], std::get<Idx>(cs)) &&
           _from_python_rec < Idx + 1 > (pys, cs);
}

template<typename... Ts>
bool from_python(py_scoped_t &py, std::tuple<Ts...> &cs)
{
    std::vector<py_scoped_t> pys;
    if (!from_python(py, pys))
        return false;

    if (pys.size() != sizeof...(Ts))
    {
        std::cerr << "expected python sequence of length " << sizeof...(Ts)
                  << ", got length " << pys.size() << std::endl;
        return false;
    }

    return _from_python_rec(pys, cs);
}

// -----

template<std::size_t Idx>
constexpr bool _to_python_rec(...)
{ return true; }

template<std::size_t Idx = 0, class ...Ts,
    class = std::enable_if_t<(Idx < sizeof...(Ts))>>
bool _to_python_rec(const std::tuple<Ts...> &cs, std::vector<py_scoped_t> &pys)
{
    return to_python(std::get<Idx>(cs), pys[Idx]) &&
           _to_python_rec<Idx + 1>(cs, pys);
}

template<typename... Ts>
bool to_python(const std::tuple<Ts...> &cs, py_scoped_t &py)
{
    std::vector<py_scoped_t> pys(sizeof...(Ts));

    if (!_to_python_rec(cs, pys))
    {
        // destruct any python references that were created prior to the error
        pys.resize(0);
    }

    return to_python(pys, py);
}

/* --------------------------------------------------------------------- */
// Helpers for enum types with an enum_map.
// Serializes to/from string.

template<typename E>
bool from_python_by_enum_map(py_scoped_t &py, E &out, E null_value,
    const char *message)
{
    std::string s;
    if (!from_python(py, s))
        return false;

    out = enum_from_str(s, null_value);
    if (out == null_value)
    {
        std::cerr << message << std::endl;
        return false;
    }
    return true;
}


} // python
} // sp2

#endif //SP2_PYTHON_CONVERSION_HPP
