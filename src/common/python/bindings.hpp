// python.hpp
//
// This module exposes an interface which does not involve any Python API types,
//  so that it does not need to include Python.h (and therefore does not force
//  consumers to worry about Python.h's crazy constraints on inclusion order)

#ifndef SP2_PYTHON_BINDINGS_HPP
#define SP2_PYTHON_BINDINGS_HPP

#include <vector>
#include <string>

namespace sp2 {
namespace python {

// Initialize the python interpreter.
//
// This should only be called once over the execution of a single program.
void initialize(const char *prog);

// Clean up the interpreter after initialize, ensuring that destructors are called, etc.
//
// This should only be called once over the execution of a single program.
int finalize();

std::vector<double> call_2d_vector_function(const char *mod_name, const char *func_name, std::vector<double> input,
        std::size_t width);

// Ensure that the given directories are in sys.path, so that the modules therein may be loaded.
//
// The paths are prepended, giving them higher priority over existing entries.
// Keep in mind that built-in modules will still take absolute top priority.
//
// Any paths already in sys.path will be moved to the front, without creating a duplicate entry.
void extend_sys_path(std::vector<std::string> dir);

}
}

#endif // SP2_PYTHON_BINDINGS_HPP
