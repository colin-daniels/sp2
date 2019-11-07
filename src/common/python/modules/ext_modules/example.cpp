#include "Python.h"

#include "common/python/modules/ext_modules.hpp"

#include "common/python/types/py_ref_t.hpp"
#include "common/python/types/py_object_t.hpp"
#include "common/python/types/as_ndarray.hpp"
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
    (void)self; // pretend to use

    // NOTE: code before wrap_cxx_logic must not throw. Stick to barebones C.
    long first;
    long second;
    const char* kw_arg_names[] = {"first", "second", nullptr};
    auto parse_ok = bool(PyArg_ParseTupleAndKeywords(args, kw, "ll:add",
        const_cast<char **>(kw_arg_names),
        &first, &second));

    if (!parse_ok)
        return nullptr; // PyErr is already set; let it propagate back up

    return wrap_cxx_logic([&] {
        long sum = first + second;
        return py_from(sum).inner().dup();
    });
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
