#ifndef SP2_PYTHON_EXT_MODULES_HPP
#define SP2_PYTHON_EXT_MODULES_HPP

// ext_modules -- C/C++ extension modules for Python

// NOTE: This file contains lists of all extension modules which must be
//       manually maintained.  Scroll down to the bottom to see it.

#include "common/python/util.hpp"

#include <functional>
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

namespace {
namespace method_def {

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

} // namespace method_def

/// Wraps the bulk of C++ logic in module functions exposed to python.
/// See the 'example' ext_module for intended usage.
///
/// The wrapped function should:
///  - Return a non-null py_scoped_t. (this is checked)
///  - Throw a std::string or const char* literal on failure.
///    This will create a Python exception.
///  - Be aware that anything else which is thrown and uncaught will also
///    converted to Python errors (though perhaps with less nice messages)
///
/// and the output of this function should be returned immediately and
///  unconditionally from the callback (since PyErr may have been set).
PyObject* wrap_cxx_logic(std::function<py_scoped_t()> func)
{
    try
    {
        py_scoped_t out = func();
        if (!out)
        {
            PyErr_SetString(PyExc_AssertionError, "wrap_cxx_logic: nullptr");
            return nullptr;
        }

        // give away ownership of the reference to the interpreter
        return out.steal();
    }
    catch (const char *e)
    {
        // This is to let you throw custom error messages directed at the python
        //  user, directly from the try block.
        PyErr_SetString(PyExc_RuntimeError, e);
        return nullptr;
    }
    catch (const std::string &e)
    {
        // This is to let you throw custom error messages directed at the python
        //  user, directly from the try block.
        PyErr_SetString(PyExc_RuntimeError, e.c_str());
        return nullptr;
    }
    catch (const std::exception &e)
    {
        std::string s =
            "An unhandled exception occurred in extension module code: ";
        s += e.what();
        PyErr_SetString(PyExc_RuntimeError, s.c_str());
        return nullptr;
    }
    catch (...)
    {
        // Getting a bit creative with 'throw', huh?
        //
        // But seriously: We need to catch EVERYTHING, no exceptions. (uh,
        // no pun intended).  Otherwise, since we are embedding the interpreter,
        // we run the risk of control flow resuming somewhere outside the Python
        // API call, thereby leaving the interpreter in an inconsistent state.
        //
        // ...I think.
        PyErr_SetString(PyExc_RuntimeError,
            ("Something unusual was thrown in extension module code."
                " That's all we know.")
        );
        return nullptr;
    }
}

} // anonymous namespace

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
