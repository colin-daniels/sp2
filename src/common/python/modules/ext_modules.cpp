#include <Python.h>

#include "ext_modules.hpp"
#include "common/python/types/py_ref_t.hpp"
#include "common/python/error.hpp"

namespace sp2 {
namespace python {
namespace ext_modules {

using namespace std;

vector<string> split(const string& s, const string& sep)
{
    vector<string> out;
    size_t start = 0;
    size_t i;
    while ((i = s.find(sep, start)) != std::string::npos) {
        out.emplace_back(s, start, i - start);
        start = i + sep.size();
    }
    out.emplace_back(s, start, s.size() - start);
    return out;
}

string join(vector<string> components, const string& sep)
{
    if (components.empty())
        return ""s;

    string last = components.back();
    components.pop_back();

    string s;
    for (const auto& c : components)
        s += c + sep;
    s += last;
    return s;
}

// Returns a (probably) unique, but undotted, name to secretly import an
// extension module under as a builtin.
std::string name_as_builtin(const char* qualified_name)
{
    auto components = split("_"s + qualified_name, ".");
    return join(components, "_");
}

enum class init_state_t {
    READY_FOR_PRE = 0,
    READY_FOR_POST = 1,
    DONE = 2
};
static auto init_state = init_state_t::READY_FOR_PRE;

void pre_py_initialize()
{
    if (Py_IsInitialized())
        throw logic_error("pre_py_initialize() called after Py_Initialize");
    if (init_state != init_state_t::READY_FOR_PRE)
        throw logic_error("pre_py_initialize() called twice");

    // horrible hack to leak temporary std::strings into c strings that will
    // outlive the python interpreter (among other things).
    //
    // this could be avoided using RAII but bleh, who cares about a couple of
    // short fixed strings produced on program initialization?
    auto leak = [](const std::string &s) -> char* {
        char* leaked = static_cast<char*>(calloc(s.size() + 1, sizeof(char)));
        std::copy(s.begin(), s.end(), leaked);
        return leaked;
    };

    // Initialize modules at the C++ level, and add them to the table of
    //  builtin modules under mangled names.
    for (auto p : all)
    {
        char* mangled_name = leak(name_as_builtin(p->qualified_name));

        if (PyImport_AppendInittab(mangled_name, p->py_initialize))
            throw runtime_error("unknown error appending to Python _inittab");
    }

    init_state = init_state_t::READY_FOR_POST;
}

/// Try to import a module through the typical machinery, or implicitly create
///  an empty one.
py_ref_t import_or_create(const char* qual_name)
{
    auto module = scope(PyImport_ImportModule(qual_name));
    if (module)
    {
        // Could be any of the following:
        //  - a legitimate, installed python module with this qualified name
        //  - something we hacked into sys.modules
        //  - an existing stub from a previous call to this function
        //
        // In any case, return it.
        return module;
    }

    // An error occurred, but what was it?
    if (!PyErr_ExceptionMatches(PyExc_ModuleNotFoundError))
        throw_on_py_err("exception importing module that appears to exist");

    // No such module; create a stub.
    PyErr_Clear();
    module = scope_dup(PyImport_AddModule(qual_name));
    throw_on_py_err("error making empty module");
    return module;
}

void post_py_initialize()
{
    if (!Py_IsInitialized())
        throw logic_error("post_py_initialize() called before Py_Initialize");
    if (init_state != init_state_t::READY_FOR_POST)
        throw logic_error("incorrect state for post_py_initialize()");

    for (auto p : all)
    {
        // During Py_Initialize, our modules will have been initialized and
        // loaded under flat, mangled names. Now make it look like they actually
        // exist at the qualified paths we wanted.
        //
        // It is not yet clear what the intended solution is for performing this
        // remapping, but the following steps make for what appears to be an
        // impenetrable illusion:

        // 1. add an entry for the qualified name to sys.modules
        string mangled_name = name_as_builtin(p->qualified_name);

        auto py_module = scope(PyImport_ImportModule(mangled_name.c_str()));
        throw_on_py_err("error importing extension module as builtin");

        auto sys_modules = scope_dup(PyImport_GetModuleDict());
        throw_on_py_err("error getting sys.modules");

        PyDict_SetItemString(sys_modules.raw(), p->qualified_name, py_module.raw());
        throw_on_py_err("error adding qualified name to sys.modules");

        // 2. ensure that each ancestor to the module can be imported,
        //    creating empty modules if need be  (import_or_create())
        // 3. equip each ancestor with an attribute for its direct child.
        auto components = split(p->qualified_name, ".");

        auto parent_path = components[0];
        auto parent = import_or_create(parent_path.c_str());
        components.erase(components.begin());

        for (auto child_name : components)
        {
            auto child_path = parent_path + "." + child_name;
            auto child = import_or_create(child_path.c_str());

            setattr(parent, child_name.c_str(), child);

            parent.destroy(); // guaranteed safe because refcount > 1
            parent = child;
            parent_path = child_path;
        }
    }

    // And that's all, folks.
    init_state = init_state_t::DONE;
};

} // namespace ext_modules
} // namespace python
} // namespace sp2
