// python.hpp
//
// This module exposes an interface which does not involve any Python API types,
//  so that it does not need to include Python.h (and therefore does not force
//  consumers to worry about Python.h's crazy constraints on inclusion order)

#ifndef SP2_PYTHON_BINDINGS_HPP
#define SP2_PYTHON_BINDINGS_HPP

#include <vector>
#include <string>

#include "phonopy/structural_mutation.hpp"

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

// Ensure that the given directories are in sys.path, so that the modules therein may be loaded.
//
// The paths are prepended, giving them higher priority over existing entries.
// Keep in mind that built-in modules will still take absolute top priority.
//
// Any paths already in sys.path will be moved to the front, without creating a duplicate entry.
void extend_sys_path(std::vector<std::string> dir);

// A highly specialized function that is inexorably tied to business logic in run_phonopy.
//
// It calls a named function in a named module (which must be available on sys.path)
// with, uh... some specific arguments, in a specific manner, and uh.... transforms the
// result back into c++ data in a specific way.
//
// Any further detail is subject to change.
sp2::structural_mutation_t call_run_phonopy_mutation_function(
    const char *mod_name, const char *func_name,
    std::vector<double> input, double lattice[3][3],
    std::vector<size_t> sc_to_prim);

}
}

#endif // SP2_PYTHON_BINDINGS_HPP
