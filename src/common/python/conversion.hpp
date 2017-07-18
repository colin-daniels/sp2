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

#include <cstdint>
#include <vector>

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
bool to_python(py_scoped_t &c, py_scoped_t &py);

bool from_python(py_scoped_t &py, long &c);
bool from_python(py_scoped_t &py, long long &c);
bool from_python(py_scoped_t &py, unsigned long &c);
bool from_python(py_scoped_t &py, unsigned long long &c);
bool from_python(py_scoped_t &py, double &c);
bool from_python(py_scoped_t &py, bool &c);
bool from_python(py_scoped_t &py, std::string &c);
bool from_python(py_scoped_t &py, std::nullptr_t &c);
bool from_python(py_scoped_t &py, py_scoped_t &c);

/* --------------------------------------------------------------------- */
// Generic conversions of vector <-> list

template <typename T>
bool to_python(const std::vector<T> & vec, py_scoped_t &list) {
    if (list) {
        throw std::logic_error("attempted to serialize into occupied py_scoped_t");
    }

    // of course, pre-allocating a list of the right size would be faster,
    // but 'new' and 'append' keeps the object in a consistent state
    list = scope(PyList_New(0));
    throw_on_py_err("unexpected error constructing python list");

    for (auto & x : vec) {
        py_scoped_t px;
        if (!to_python(x, px)) {
            // ...and here's why it's nice to keep 'list' in a consistent state
            list.destroy();
            return false;
        }

        PyList_Append(list.raw(), px.raw());
        if (print_on_py_err()) {
            list.destroy();
            throw std::runtime_error("unexpected error appending to python list");
        }
    }

    return true;
}

template <typename T>
bool from_python(py_scoped_t &o, std::vector<T> & vec) {
    auto dummy = std::move(vec); // evict any existing contents

    // NOTICE
    // Although the docs for PySequence_Fast say that it takes a sequence,
    // it actually accepts iterables as well.
    //
    // It's basically:  x if isinstance(x, (list, tuple)) else list(x)
    //
    // Both of these cases are fine; but we'll check for them explicitly to
    // ensure that any errors from PySequence_Fast are not related to them.
    bool is_iterable = bool(scope(PyObject_GetIter(o.raw())));
    bool is_sequence = PySequence_Check(o.raw());
    bool is_iterator = PyIter_Check(o.raw());

    if (!(is_iterable || is_sequence)) {
        return false;
    }

    // Also, we'd rather NOT accept iterators, which would get consumed the
    // first time we fail to deserialize them and then fail to deserialize
    // as anything else.
    if (is_iterator) {
        return false;
    }

    // With that out of the way, any errors from PySequence_Fast are
    // legitimate concerns.
    py_scoped_t seq = scope(PySequence_Fast(o.raw(), "expected a... bah, you'll never see this."));
    throw_on_py_err();

    Py_ssize_t size = PySequence_Fast_GET_SIZE(o.raw());
    vec.reserve(size);
    for (Py_ssize_t i=0; i < size; i++) {
        auto py_item = scope(PySequence_Fast_GET_ITEM(o.raw(), i));

        T item;
        if (!from_python(py_item.dup(), item)) {
            return false;
        }

        vec.push_back(item);
    }

    return true;
}

/* --------------------------------------------------------------------- */
// Generic conversions of ndarray_serialize_t <-> ndarray

/// Helper type to obtain the numpy dtype integral constant associated
/// with a data type.
template <typename T>
struct numpy_dtype {};

#define SPECIALIZE_NUMPY_DTYPE(T,VAL) \
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
template <typename T, int DTYPE = numpy_dtype<T>::value>
bool to_python(const ndarray_serialize_t<T> &c, py_scoped_t &py)
{

    // copy data into a brand new array object.
    std::vector<npy_intp> shape(c.shape().begin(), c.shape().end());
    auto arr = scope(PyArray_SimpleNew(c.ndim(), shape.data(), DTYPE));
    if (print_on_py_err())
        return false;

    auto arr_data = (double *)PyArray_DATA((PyArrayObject *)arr.raw());
    throw_on_py_err("error accessing numpy array data");
    copy(c.data().begin(), c.data().end(), arr_data);

    py = std::move(arr);
    return true;
}

// NOTE: DTYPE doubles as SFINAE
template <typename T, int DTYPE = numpy_dtype<T>::value>
bool from_python(py_scoped_t &py, ndarray_serialize_t<T> &c)
{

    // Force the array into a contiguous layout if it isn't.
    int min_depth = 0; // ignore
    int max_depth = 0; // ignore
    auto contiguous = scope(PyArray_ContiguousFromAny(py.raw(), DTYPE, min_depth, max_depth));
    if (print_on_py_err()) {
        return false;
    }

    auto arr = (PyArrayObject *)contiguous.raw();
    // NOTE: none of these set the python error state
    auto arr_data = (double *)PyArray_DATA(arr);
    size_t arr_size = PyArray_SIZE(arr);
    size_t arr_ndim = PyArray_NDIM(arr);
    auto* arr_dims = PyArray_DIMS(arr);

    std::vector<double> data(arr_data, arr_data + arr_size);
    std::vector<size_t> shape(arr_dims, arr_dims + arr_ndim);
    c = ndarray_serialize_t<double>(data, shape);
    return true;
}

} // python
} // sp2

#endif //SP2_PYTHON_CONVERSION_HPP
