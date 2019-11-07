#ifndef SP2_PYTHON_ENVIRONMENT_HPP
#define SP2_PYTHON_ENVIRONMENT_HPP

// environment.hpp:
//
//   Initializing and configuring the python interpreter.

#include <vector>
#include <string>

namespace sp2 {
namespace python {

/// Handles state of the global python interpreter through RAII.
class environment
{
    wchar_t *py_allocated_program = nullptr;

    void ensure_run_once();
    int finalize();

public:
    /// Initialize the python interpreter.  Expects to receive argv[0].
    ///
    /// Only one environment may be constructed over the course of a program's
    /// execution.
    environment(const char *prog);

    /// Clean up the interpreter, ensuring that __del__ methods are called, etc.
    ~environment();
};

/// Ensure that the given directories are in sys.path,
/// so that the modules therein may be loaded.
///
/// The paths are prepended, giving them higher priority over existing entries.
/// Keep in mind that built-in modules will still take absolute top priority.
/// Any paths already in sys.path will be moved to the front, without creating
/// a duplicate entry.
void extend_sys_path(std::vector<std::string> dir);

} // namespace python
} // namespace sp2

#endif // SP2_PYTHON_ENVIRONMENT_HPP
