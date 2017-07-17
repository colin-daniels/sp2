#include "Python.h" // must be first include

#include "conversion.hpp"

namespace sp2 {
namespace python {

void refuse_to_evict(py_scoped_t &p) {
    if (p) {
        throw std::logic_error("attempted to serialize into already-occupied py_scoped_t");
    }
}

// basic implementation given a function with these properties:
//  - type (Args...) -> PyObject*
//  - Returns a new reference on success
//  - Returns NULL and sets the python error state on failure
template <typename F, typename... Args>
bool to_python_simple_conv(py_scoped_t &py, F func, Args... args) {
    refuse_to_evict(py);

    py = scope(func(args...));
    return !print_on_py_err();
};

bool to_python(const nullptr_t &c, py_scoped_t &py) {
    refuse_to_evict(py);

    py = scope_dup(Py_None);
    return true;
}

bool to_python(const long &c, py_scoped_t &py) {
    return to_python_simple_conv(py, PyLong_FromLong, c);
}

bool to_python(const unsigned long &c, py_scoped_t &py) {
    return to_python_simple_conv(py, PyLong_FromUnsignedLong, c);
}

bool to_python(const long long &c, py_scoped_t &py) {
    return to_python_simple_conv(py, PyLong_FromLongLong, c);
}

bool to_python(const unsigned long long &c, py_scoped_t &py) {
    return to_python_simple_conv(py, PyLong_FromUnsignedLongLong, c);
}

bool to_python(const double &c, py_scoped_t &py) {
    return to_python_simple_conv(py, PyFloat_FromDouble, c);
}

bool to_python(const bool &c, py_scoped_t &py) {
    refuse_to_evict(py);

    if (c) {
        py = scope_dup(Py_True);
    } else {
        py = scope_dup(Py_False);
    }
    return true;
}

bool to_python(const std::string &utf8, py_scoped_t &py) {
    return to_python_simple_conv(py,
            PyUnicode_DecodeUTF8, utf8.c_str(), utf8.size(), "strict");
}

bool to_python(py_scoped_t &c, py_scoped_t &py) {
    refuse_to_evict(py);

    py = c.dup();
}

/* --------------------------------------------------------------------- */

// FIXME these implementations do not differentiate between unrecoverable
//       failures and tolerable failures (e.g. wrong type)

bool from_python(py_scoped_t py, long &c) {
    c = PyLong_AsLong(py.raw());
    return !print_on_py_err();
}

bool from_python(py_scoped_t py, unsigned long &c) {
    c = PyLong_AsUnsignedLong(py.raw());
    return !print_on_py_err();
}

bool from_python(py_scoped_t py, long long &c) {
    c = PyLong_AsLongLong(py.raw());
    return !print_on_py_err();
}

bool from_python(py_scoped_t py, unsigned long long &c) {
    c = PyLong_AsUnsignedLongLong(py.raw());
    return !print_on_py_err();
}

bool from_python(py_scoped_t py, double &c) {
    c = PyFloat_AsDouble(py.raw());
    return !print_on_py_err();
}

bool from_python(py_scoped_t py, bool &c) {
    if (py.raw() == Py_True) {
        c = true;
        return true;

    } else if (py.raw() == Py_False) {
        c = false;
        return true;
    }

    return false;
}

bool from_python(py_scoped_t py, nullptr_t &c) {
    c = nullptr; // *shrug*
    return py.raw() == Py_None;
}

bool from_python(py_scoped_t py, std::string &c) {
    std::string dummy = move(c); // evict contents

    // From the Python API docs:
    //     This caches the UTF-8 representation of the string in the
    //     Unicode object, and subsequent calls will return a pointer
    //     to the same buffer.
    //     The caller is not responsible for deallocating the buffer.
    char* buf = PyUnicode_AsUTF8(py.raw());
    if (print_on_py_err()) {
        return false;
    }
    c = std::string(buf);
    return true;
}

bool from_python(py_scoped_t py, py_scoped_t &c) {
    refuse_to_evict(c);

    c = py.dup();
}


} // namespace python
} // namespace sp2
