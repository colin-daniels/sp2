#include <Python.h> // must be first include

#include "fake_modules.hpp"

#include "error.hpp"
#include "py_scoped_t.hpp"

#include <string>

using std::string;

void sp2::python::fake_modules::initialize()
{
    for (auto p : fake_modules::all)
    {
        // This just shows up in debugging info... probably.
        auto fake_path = string() + "sp2/common/python/fake_modules/" + p->meta.name + ".cpp";

        // This identifies the module in sys.modules
        auto fake_qualified_name = string() + fake_modules::PACKAGE + "." + p->meta.name;

        py_scoped_t code = scope(Py_CompileString(p->meta.text, fake_path.c_str(), Py_file_input));
        throw_on_py_err();
        py_scoped_t module = scope(PyImport_ExecCodeModule(fake_qualified_name.c_str(), code.raw()));
        throw_on_py_err();

        p->module = std::move(module);
    }
}

void sp2::python::fake_modules::finalize()
{
    for (auto p : fake_modules::all)
        p->module.destroy();
}
