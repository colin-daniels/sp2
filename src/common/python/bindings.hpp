// python.hpp
//
// This module exposes an interface which does not involve any Python API types,
//  so that it does not need to include Python.h (and therefore does not force
//  consumers to worry about Python.h's crazy constraints on inclusion order)

#ifndef SP2_PYTHON_BINDINGS_HPP
#define SP2_PYTHON_BINDINGS_HPP

#include <vector>
#include <string>
#include <memory>

#include "phonopy/structural_mutation.hpp"

namespace sp2 {
namespace python {

/// An opaque, pImpl-style wrapper around a python pointer.
struct py_opaque_t
{
    struct impl_t;
    std::unique_ptr<impl_t> impl;

    py_opaque_t(std::unique_ptr<impl_t>&&);

    // the default implementations for these suffice, but we need to make sure
    // they are generated in a location where 'impl_t' is fully defined.
    py_opaque_t();
    py_opaque_t(py_opaque_t&&);
    py_opaque_t& operator=(py_opaque_t&&);
    ~py_opaque_t();
    py_opaque_t(const py_opaque_t&) = delete;
    py_opaque_t& operator=(const py_opaque_t&) = delete;
};

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

/// Obtain module.attr, which must exist (else an exception is thrown)
///
/// This implicitly imports the module through the standard python
/// module-loading machinery, which caches imports.
///
/// In other words, this MAY invoke module.__init__() (and therefore one should
///  be prepared for arbitrary python errors), but it will only do so once.
py_opaque_t lookup_required(const char *mod_name, const char *attr);

/// Obtain module.attr, or a NULL pointer.
///
/// Like 'lookup_required', this may or may not invoke 'module.__init__()'.
py_opaque_t lookup_optional(const char *mod_name, const char *attr);

/// Highly specialized functions that are inexorably tied to business
/// logic in run_phonopy.
namespace run_phonopy {

// It calls a named function in a named module (which must be available on sys.path)
// with, uh... some specific arguments, in a specific manner, and uh.... transforms the
// result back into c++ data in a specific way.
//
// Any further detail is subject to change.
sp2::python::py_opaque_t make_param_pack(
    std::vector<double> carts, const double lattice[3][3],
    std::vector<size_t> sc_to_prim);

// function must NOT be NULL
sp2::structural_mutation_t call_mutate(
    const py_opaque_t &function,
    const py_opaque_t &param_pack);

// function can be NULL
sp2::structural_mutation_t call_apply(
    const py_opaque_t &function,
    const py_opaque_t &mutation,
    const py_opaque_t &param_pack);

// function can be NULL
void call_accept(
    const py_opaque_t &function,
    const py_opaque_t &mutation,
    const py_opaque_t &param_pack);

// function must NOT be NULL
py_opaque_t call_generate(
    const py_opaque_t &function,
    const py_opaque_t &param_pack);

} // namespace run_phonopy
} // namespace python
} // namespace sp2

#endif // SP2_PYTHON_BINDINGS_HPP
