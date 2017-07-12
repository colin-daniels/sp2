
#ifndef SP2_PYTHON_BINDINGS_HPP
#define SP2_PYTHON_BINDINGS_HPP

// Must be included before standard headers
#include <Python.h>

#include <stdexcept>
#include <vector>
#include <object.h>

/* --------------------------------------------------------------------- */

namespace sp2 {
namespace python {

// A scoped reference to a python object that uses RAII to handle Py_DECREF.
// This makes it somewhat easier to reason about exception safety.
//
// The contained object may be NULL.

class py_scoped_t
{
    PyObject *obj;

    // null constructor
    py_scoped_t() {};

    // private constructor
    explicit py_scoped_t(PyObject * o);

    // explicit constructor from a borrowed ref, which makes a new reference
    friend py_scoped_t scope_dup(PyObject * o);

    // explicit constructor from a new ref
    friend py_scoped_t scope(PyObject * o);

public:

    explicit operator bool() const;

    // No copying.  Use dup() to make the incref explicit.
    py_scoped_t& operator=(const py_scoped_t& other) = delete;
    py_scoped_t(const py_scoped_t& other) = delete;

    // Moving.
    py_scoped_t(py_scoped_t&& other);

    py_scoped_t& operator=(py_scoped_t&& other);

    ~py_scoped_t();

    // Increment the refcount and return a new scoped reference.
    py_scoped_t dup();

    // Borrow the reference without touching the refcount.
    //
    // This is the appropriate method for interfacing with most Python APIs.
    PyObject * raw() const;

    // Leak the reference, preventing the DECREF that would otherwise occur at scope exit.
    // The scoped reference will become NULL.
    //
    // Necessary for working with Python API functions that steal references,
    // such as PyTuple_SetItem.
    PyObject * steal();
};

py_scoped_t scope_dup(PyObject * o);
py_scoped_t scope(PyObject * o);

std::vector<double> call_vector_function(const char *mod_name, const char *func_name, std::vector<double> input);

py_scoped_t call_module_function_kw(const char *mod_name, const char *func_name, py_scoped_t kw);

py_scoped_t call_module_function(const char *mod_name, const char *func_name, py_scoped_t args, py_scoped_t kw);

// Initialize the python interpreter.
//
// This should only called once over the execution of a single program.
void initialize(const char *prog);

// Clean up the interpreter after initialize, ensuring that destructors are called, etc.
//
// This should only called once over the execution of a single program.
int finalize();



/* SOMEDAY
static PyObject *
positions_as_ndarray(std::vec<double> &positions)
{
    npy_intp dim[2] = {3, positions.size()/3};
    return PyArray_SimpleNewFromData(
            2, dim, // ndim, dim,
            NPY_DOUBLE, positions.data());
}
*/

// Ensure that a given directory is in sys.path, so that the modules therein may be loaded.
//
// This function actively checks for and avoids adding duplicate entries.
// A return value of 'false' means sys.path was not modified.
bool add_to_sys_path(const char* dir);

std::vector<double> flat_2d_from_py_sequence(size_t ncol, PyObject *o);

py_scoped_t py_list_from_flat_2d(size_t nrow, size_t ncol, double *data);

std::vector<double> vec_from_py_sequence(PyObject *o);

py_scoped_t py_list_from_vec(std::vector<double> v);

// get an object's repr(), mostly for debug purposes
std::wstring repr(py_scoped_t o);

// get an object's str(), mostly for debug purposes
std::wstring str(py_scoped_t o);

}
}
#endif // SP2_PYTHON_BINDINGS_HPP
