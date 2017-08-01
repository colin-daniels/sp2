#include "Python.h"

#include "common/python/ext_modules.hpp"

#include "common/python/util.hpp"
#include "common/python/numpy_util.hpp"
#include "common/python/conversion.hpp"

namespace sp2 {
namespace python {
namespace ext_modules {
namespace example {

// NOTE: Pretend this says 'sp2.example'.
//       This is mangled only because there is also a fake_module example.
const char *QUALIFIED_NAME = "sp2._example_ext";

using namespace std;

PyObject* py_add(PyObject *self, PyObject *args, PyObject *kw)
{
    long first;
    long second;
    const char* kw_arg_names[] = {"first", "second", nullptr};
    auto parse_ok = bool(PyArg_ParseTupleAndKeywords(args, kw, "ll:add",
        const_cast<char **>(kw_arg_names),
        &first, &second));

    if (!parse_ok)
        return NULL; // PyErr is already set; let it propagate back up

    // wrap everything in a catchall to avoid unwinding through Python
    // where it would bypass cleanup code.  This is especially important
    // given how we are embedding the interpreter, as otherwise, such exceptions
    // could end up getting caught.
    try
    {
        long sum = first + second;
        py_scoped_t py_sum = to_python_strict(sum);

        return py_sum.steal(); // returned references must be leaked
    }
    catch (std::string e)
    {
        // This is to let you throw custom error messages directed at the python
        //  user, directly from the try block.
        PyErr_SetString(PyExc_RuntimeError, e.c_str());
        return NULL;
    }
    catch (const exception &e)
    {
        // if all else fails, convert C++ exceptions to python exceptions
        PyErr_SetString(PyExc_RuntimeError,
            ("An unhandled exception occurred in extension module code: "s
            + e.what()).c_str()
           );
        return NULL;
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
        return NULL;
    }
}

PyMethodDef py_method_defs[] = {
    method_def::args_kw("add", py_add, "Add two integers."),
    method_def::END
};

PyModuleDef py_module_def = {
    PyModuleDef_HEAD_INIT, // PyModuleDef_Base m_base
    QUALIFIED_NAME,        // char* m_name
    NULL,                  // char* m_doc
    0,                     // Py_ssize_t m_size
    py_method_defs,        // PyMethodDef* m_methods
    NULL,                  // PyModuleDef_Slot* m_slots
    NULL,                  // traverseproc m_traverse
    NULL,                  // inquiry m_clear
    NULL                   // freefunc m_free
};

PyObject* py_initialize()
{ return PyModuleDef_Init(&py_module_def); }

ext_module_t ext_module { QUALIFIED_NAME, py_initialize };

} // namespace example
} // namespace ext_modules
} // namespace python
} // namespace sp2
