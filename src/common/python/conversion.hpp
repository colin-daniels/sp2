#ifndef SP2_PYTHON_CONVERSION_HPP
#define SP2_PYTHON_CONVERSION_HPP

#ifndef Py_PYTHON_H
#error This module is designed for internal use by source files in common/python, \
and is not safe to use unless Python.h is already included before all standard \
library headers.
#endif // Py_PYTHON_H

#include "common/python/util.hpp"

#include <vector>

namespace sp2 {
namespace python {

#warning these python conversion functions are all WILDLY untested

bool to_python(const long &c, py_scoped_t &py);
bool to_python(const long long &c, py_scoped_t &py);
bool to_python(const unsigned long &c, py_scoped_t &py);
bool to_python(const unsigned long long &c, py_scoped_t &py);
bool to_python(const double &c, py_scoped_t &py);
bool to_python(const bool &c, py_scoped_t &py);
bool to_python(const std::string &c, py_scoped_t &py);
bool to_python(const std::nullptr_t &c, py_scoped_t &py);
bool to_python(py_scoped_t &c, py_scoped_t &py);

bool from_python(py_scoped_t py, long &c);
bool from_python(py_scoped_t py, long long &c);
bool from_python(py_scoped_t py, unsigned long &c);
bool from_python(py_scoped_t py, unsigned long long &c);
bool from_python(py_scoped_t py, double &c);
bool from_python(py_scoped_t py, bool &c);
bool from_python(py_scoped_t py, std::string &c);
bool from_python(py_scoped_t py, std::nullptr_t &c);
bool from_python(py_scoped_t py, py_scoped_t &c);

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
bool from_python(py_scoped_t o, std::vector<T> & vec) {
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

} // python
} // sp2

#endif //SP2_PYTHON_CONVERSION_HPP
