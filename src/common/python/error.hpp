#ifndef SP2_PYTHON_ERROR_HPP
#define SP2_PYTHON_ERROR_HPP

namespace sp2 {
namespace python {

/// Checks for python exceptions via the PyErr API
/// and turns them into C++ exceptions.
void throw_on_py_err(const char *msg);

void throw_on_py_err();

/// Checks for a python exception and prints it to stderr, clearing the error.
///
/// Returns true if an error was found.
bool print_on_py_err();

} // namespace python
} // namespace sp2
#endif // SP2_PYTHON_ERROR_HPP
