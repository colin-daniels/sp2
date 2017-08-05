#ifndef SP2_PYTHON_ERROR_HPP
#define SP2_PYTHON_ERROR_HPP

#include "common/python/types/py_object_t_body.hpp"

#include <array>
#include <string>
#include <stdexcept>
#include <memory>

namespace sp2 {
namespace python {

static const char * DEFAULT_PY_ERROR_MSG = "An exception was thrown in Python.";

typedef std::array<py_object_t, 3> py_exc_info;

/// Represents a Python exception, possibly with a string message indicating
///  where the exception occurred in C++ code.
///
/// For best diagnostics, catch this type explicitly and use 'print_to_stderr',
/// then possibly rethow as 'std::runtime_error'. The polymorphic 'what()'
/// method will not contain information about the python error.
class py_error: public std::runtime_error
{
    std::shared_ptr<py_exc_info> info;

    // invariant: info[0] (the exception type) is not NULL.

public:

    py_error(const py_exc_info& info, const std::string &message)
        : std::runtime_error(message), info(std::make_shared<py_exc_info>(info))
    {
        const auto& arr = *(this->info);
        if (!arr[0])
            throw std::logic_error("created py_error with NULL exc_type");
    }

    /// The type of the exception.  Always present.
    const py_object_t& exc_type() const { return (*this->info)[0]; }
    /// The value of the exception.  Can be NULL.
    const py_object_t& exc_value() const { return (*this->info)[1]; }
    /// A python stack trace associated with the exception.  Can be NULL.
    const py_object_t& exc_traceback() const { return (*this->info)[2]; }

    /// Print the typical python exception info we all know and love,
    /// tracebacks included, directly to stderr.
    ///
    /// Need it to go somewhere else? Well, tough. This is what CPython exposes.
    ///
    /// Returns a reference so you can then invoke 'rethrow()' if need be.
    const py_error& print_to_stderr() const;

    /// Rethrow as 'std::runtime_error'.
    void rethrow() const;
};

/// Extracts the global python error state, clearing it.
///
/// Returns NULL if there is no error.
std::unique_ptr<py_error> grab_py_err(const char* msg);
static std::unique_ptr<py_error> grab_py_err()
{ return grab_py_err(DEFAULT_PY_ERROR_MSG); }

/// Set the global python error state.
void restore_py_err(const py_error& err);

/// Checks for python exceptions via the PyErr API
/// and turns them into 'py_error' exceptions.
void throw_on_py_err(const char *msg);
static void throw_on_py_err()
{ throw_on_py_err(DEFAULT_PY_ERROR_MSG); }

/// Checks for a python exception and prints it to stderr, clearing the error.
///
/// Returns true if an error was found.
bool print_on_py_err();

} // namespace python
} // namespace sp2
#endif // SP2_PYTHON_ERROR_HPP
