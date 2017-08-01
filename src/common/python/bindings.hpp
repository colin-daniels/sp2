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

#include "common/minimize/metropolis_enums.hpp"

namespace sp2 {
namespace python {

struct py_scoped_t;

/// An opaque, pImpl-style wrapper around a python pointer
/// for safe encapsulation even in code that does not import
/// the cPython headers.
class py_opaque_t
{
public:
    /// pImpl target type.
    typedef py_scoped_t impl_t;

private:
    std::shared_ptr<impl_t> _impl;

public:

    py_scoped_t& inner();
    const py_scoped_t& inner() const;

    py_opaque_t(std::shared_ptr<py_scoped_t>&&);

    py_opaque_t();
    py_opaque_t(py_opaque_t&&);
    py_opaque_t& operator=(py_opaque_t&&);
    ~py_opaque_t();
    py_opaque_t(const py_opaque_t&);
    py_opaque_t& operator=(const py_opaque_t&);

private:
    void debug_null() const;
};

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
///
/// Any paths already in sys.path will be moved to the front,
/// without creating a duplicate entry.
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

/// Conflict-resolution strategies for merge_dictionaries.
enum class merge_strategy : int
{
    /// Resolve key conflicts by taking the first dict's value.
    USE_FIRST = 0,
    /// Resolve key conflicts by taking the second dict's value.
    USE_SECOND = 1,
    /// Don't resolve key conflicts; throw a runtime_exception
    ERROR = 2
};

/// Perform a union-like operation on two python dictionaries that produces a
/// dict with all of their (key, value) pairs.
py_opaque_t merge_dictionaries(
    const py_opaque_t &a, const py_opaque_t &b,
    merge_strategy strategy = merge_strategy::USE_SECOND);

/*----------------------------------------------------------------------------*/

/// Highly specialized functions that are inexorably tied to business
/// logic in run_phonopy.
namespace run_phonopy
{

/// Produce extra structural_metropolis kw args for run_phonopy.
///
/// These will not conflict with any of structural_metropolis' own kw args.
py_opaque_t make_extra_kw(std::vector<size_t> sc_to_prim);

} // namespace run_phonopy

/*----------------------------------------------------------------------------*/

/// Functions that interface with the Python interpreter to run
/// callbacks for structural_metropolis.
namespace structural_metropolis
{

/// Produce the standard kw arguments that are given to the majority of
/// structural_metropolis' python callbacks.
py_opaque_t make_param_pack(std::vector<double> carts,
    const double lattice[3][3], std::vector<double> force);

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

// function can be NULL
bool call_is_repeatable(
    const py_opaque_t &function,
    const py_opaque_t &mutation,
    const py_opaque_t &param_pack);

// function can be NULL
py_opaque_t call_scale(
    const py_opaque_t &function,
    const py_opaque_t &mutation,
    double factor,
    const py_opaque_t &param_pack);

} // namespace structural_metropolis

/*----------------------------------------------------------------------------*/

} // namespace python
} // namespace sp2

#endif // SP2_PYTHON_BINDINGS_HPP
