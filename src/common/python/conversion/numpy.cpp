#include "Python.h" // must be first include

#define SP2_NDARRAY_INSTANTIATIONS_VISIBLE
#include "common/python/conversion/numpy.hpp"

#include "common/python/include_numpy.hpp"
#include "common/python/types/py_scoped_t_body.hpp"

namespace sp2 {
namespace python {

/// Helper type to obtain the numpy dtype integral constant associated
/// with a data type.
template<typename T>
struct numpy_dtype {};

#define SP2_SPECIALIZE_NUMPY_DTYPE(T, VAL) \
template<> struct numpy_dtype<T> : std::integral_constant<int, VAL> {}

// NOTE: Be sure to keep this and SP2_FOR_EACH_NUMPY_DTYPE in sync
SP2_SPECIALIZE_NUMPY_DTYPE(bool, NPY_BOOL);
SP2_SPECIALIZE_NUMPY_DTYPE(int8_t, NPY_INT8);
SP2_SPECIALIZE_NUMPY_DTYPE(int16_t, NPY_INT16);
SP2_SPECIALIZE_NUMPY_DTYPE(int32_t, NPY_INT32);
SP2_SPECIALIZE_NUMPY_DTYPE(int64_t, NPY_INT64);
SP2_SPECIALIZE_NUMPY_DTYPE(uint8_t, NPY_UINT8);
SP2_SPECIALIZE_NUMPY_DTYPE(uint16_t, NPY_UINT16);
SP2_SPECIALIZE_NUMPY_DTYPE(uint32_t, NPY_UINT32);
SP2_SPECIALIZE_NUMPY_DTYPE(uint64_t, NPY_UINT64);
SP2_SPECIALIZE_NUMPY_DTYPE(float, NPY_FLOAT);
SP2_SPECIALIZE_NUMPY_DTYPE(double, NPY_DOUBLE);

template<typename T>
bool to_python_by_ndarray(const as_ndarray_t<T> &c, py_scoped_t &py)
{
    int dtype = numpy_dtype<T>::value;
    // copy data into a brand new array object.
    std::vector<npy_intp> shape(c.shape().begin(), c.shape().end());
    auto arr = scope(PyArray_SimpleNew(c.ndim(), shape.data(), dtype));
    if (print_on_py_err())
        return false;

    auto arr_data = (T *) PyArray_DATA((PyArrayObject *) arr.raw());
    throw_on_py_err("error accessing numpy array data");
    std::copy(c.data().begin(), c.data().end(), arr_data);

    py = std::move(arr);
    return true;
}

template<typename T>
bool from_python_by_ndarray(const py_scoped_t &py, as_ndarray_t<T> &c)
{
    int dtype = numpy_dtype<T>::value;
    // Force the array into a contiguous layout if it isn't.
    auto contiguous = scope(PyArray_ContiguousFromAny(py.raw(), dtype, 0, 0));
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

/* -------------------------------------------------------------------------- */
// Generate explicit instantiations of the template.

#define MAC(T) \
template bool to_python_by_ndarray<T>(const as_ndarray_t<T> &c, py_scoped_t &py); \
template bool from_python_by_ndarray<T>(const py_scoped_t &py, as_ndarray_t<T> &c)

SP2_FOR_EACH_NUMPY_DTYPE(MAC);
#undef MAC

} // namespace python
} // namespace sp2
