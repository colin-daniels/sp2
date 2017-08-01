#ifndef SP2_PYTHON_EXT_MODULES_HPP
#define SP2_PYTHON_EXT_MODULES_HPP

// ext_modules -- C/C++ extension modules for Python

// NOTE: This file contains lists of all extension modules which must be
//       manually maintained.  Scroll down to the bottom to see it.

#include <vector>
#include <tuple>

namespace sp2 {
namespace python {

/// Defines a Python extension module, i.e. a module available for import
///  in Python which is implemented in C/C++.
struct ext_module_t
{
    /// The package-qualified name for the module, i.e. what one would write
    /// in an 'import' statement in python.
    const char *qualified_name;

    /// This is given directly to PyImport_AppendInittab as a function pointer.
    /// Extension modules are free to choose between single-phase initialization
    /// (PyModuleDef_Init) and multi-phase (PyModule_Create)
    PyObject* (*py_initialize)();
};

} // namespace python
} // namespace sp2

/*----------------------------------------------------------------------------*/

namespace sp2 {
namespace python {
namespace ext_modules {

/*----------------------------------------------------------------------------*/
// outward-facing API

/// Performs some initialization BEFORE Py_Initialize.
void pre_py_initialize();

/// Performs some more initialization AFTER Py_Initialize.
void post_py_initialize();

/*----------------------------------------------------------------------------*/
// internal utils, defined here for convenience

namespace method_def {
namespace {

/// Helper method to make a PyMethodDef for a function that takes keyword args.
constexpr PyMethodDef
args_kw(const char *name, PyCFunctionWithKeywords func, const char *doc)
{
    // NOTE: The reinterpret_cast is from a 3-ary function pointer to a 2-ary
    //       function pointer, because that's just how PyMethodDef was defined.
    //       You are, in fact, expected to do this.
    //       https://docs.python.org/3/c-api/structures.html#c.PyMethodDef
    return {
        name, reinterpret_cast<PyCFunction>(func),
        METH_VARARGS | METH_KEYWORDS, doc
    };
}

/// Define a method that just takes *args.
constexpr PyMethodDef args(const char *name, PyCFunction func, const char *doc)
{ return {name, func, METH_VARARGS, doc}; }

/// The sentinel entry that must appear at the end of a PyMethodDef list.
constexpr PyMethodDef END = {nullptr, nullptr, 0, nullptr};

} // anonymous namespace
} // namespace py_method_def
/*----------------------------------------------------------------------------*/
// module list
//
// Please keep the list of namespaces and the vector 'all' in sync.

namespace example { extern ext_module_t ext_module; }

static std::vector<ext_module_t*> all {
    &example::ext_module
};

/*----------------------------------------------------------------------------*/
} // namespace ext_modules
} // namespace python
} // namespace sp2

#endif // SP2_PYTHON_EXT_MODULES_HPP
