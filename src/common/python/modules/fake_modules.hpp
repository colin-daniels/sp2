#ifndef SP2_PYTHON_FAKE_MODULES_HPP
#define SP2_PYTHON_FAKE_MODULES_HPP

// Fake python modules embedded in static strings.
//
// These have the advantage of not having to worry about python package
//  distribution and versioning issues.
//
// They have the disadvantage of EVERYTHING ELSE ABOUT THEM,
// but also namely:
//  * Maintaining scripts embedded in giant strings
//  * From the documentation, it doesn't look like these get cached
//    in 'sys' quite the same way that modules imported by the normal
//    mechanisms do. (thus we have a function to import them all once
//    during init)

#include "common/python/util.hpp"

#include <vector>

namespace sp2 {
namespace python {

struct fake_module_template_t {
    // Name of the module, which will become the part after the dot
    //  in its package-qualified name in Python.
    const char * name;
    // Textual content of a python script representing the
    // module's __init__.
    const char * text;
};

// A loaded fake_module.
//
// Instances of this are declared globally because it is not yet
// clear whether modules defined using PyImport_ExecCodeModule could
// be safely loaded more than once, anyways.
struct fake_module_t {
    fake_module_template_t meta;

    // Python module object.
    // This is NULL until initialize_fake_modules(),
    // and after finalize_fake_modules().
    py_scoped_t module;
};

namespace fake_modules {
    static const char * PACKAGE = "_sp2_fake";

    // These lists must be maintained in parallel:

    // 1) the module templates, each defined in their own file.
    extern fake_module_template_t mutation_helper_template;

    // 2) objects representing the loaded modules.
    //    Definitions are alongside the template, using NULL py_scoped_t() for the module.
    extern fake_module_t mutation_helper;

    // 3) A list of all modules so we can more easily iterate over them elsewhere.
    static const std::vector<fake_module_t*> all = {
        &mutation_helper
    };

} // namespace fake_modules

void initialize_fake_modules();
void finalize_fake_modules();

} // namespace python
} // namespace sp2

#endif //SP2_PYTHON_FAKE_MODULES_HPP
