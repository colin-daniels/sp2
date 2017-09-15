#include <Python.h> // must be first include

#include "fake_modules.hpp"

#include "common/python/error.hpp"
#include "common/python/types/py_ref_t.hpp"

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

        py_ref_t code = scope(Py_CompileString(p->meta.text, fake_path.c_str(), Py_file_input));
        throw_on_py_err();
        py_ref_t module = scope(PyImport_ExecCodeModule(fake_qualified_name.c_str(), code.raw()));
        throw_on_py_err();

        p->module = std::move(opaque(module));
    }
}

void sp2::python::fake_modules::finalize()
{
    for (auto p : fake_modules::all)
        p->module.destroy();
}

sp2::python::fake_module_t::~fake_module_t()
{
    // FIXME EVIL HORRIBLE HACK
    //
    // Deliberately leak the reference.
    // Fake modules are all allocated statically, and this makes things
    //  get ugly if C code invokes exit(), since they get destructed after
    //  the interpreter has been shut down.
    //
    // The proper solution is to find a way to scope their lifetime
    //  so that it isn't static, and ends before the interpreter closes.
    this->module.inner().steal();
}
