#ifndef SP2_PYTHON_FAKE_MODULES_HPP
#define SP2_PYTHON_FAKE_MODULES_HPP

// fake_modules -- Pure python modules, embedded in static strings.

// NOTE: This file contains lists of all fake modules which must be
//       manually maintained.  Scroll down to the bottom to see it.

// Python modules embedded as text in the binary have the advantage of not
// having to worry about python package distribution and versioning issues,
// but also have the disadvantage of EVERYTHING ELSE ABOUT THEM.
//
// Namely:
//  * Maintaining scripts embedded in giant strings
//  * From the documentation, it doesn't look like these get cached
//    in 'sys' quite the same way that modules imported by the normal
//    mechanisms do. (thus we have a function to import them all once
//    during init)
//
// Nonetheless, it seems workable for now.

#include "common/python/types/py_scoped_t_body.hpp"

#include <vector>

/*----------------------------------------------------------------------------*/
// data types

namespace sp2 {
namespace python {

struct fake_module_template_t
{
    /// Name of the module, which will become the part after the dot
    ///  in its package-qualified name in Python.
    const char * name;
    /// Textual content of a python script representing the
    /// module's __init__.
    const char * text;
};

/// A loaded fake_module.
///
/// Instances of this are declared globally because it is not yet
/// clear whether modules defined using PyImport_ExecCodeModule could
/// be safely loaded more than once, anyways.
///
/// In this manner, fake modules are effectively singletons.
struct fake_module_t
{
    fake_module_template_t meta;

    /// Python module object;
    ///
    /// This is NULL until initialize_fake_modules(),
    /// and after finalize_fake_modules().
    py_scoped_t module;

    fake_module_t(fake_module_template_t meta)
            : meta(meta)
    {}
};

} // namespace python
} // namespace sp2

/*----------------------------------------------------------------------------*/
// Outward-facing API.

namespace sp2 {
namespace python {
namespace fake_modules {

void initialize();
void finalize();

} // namespace fake_modules
} // namespace python
} // namespace sp2

/*----------------------------------------------------------------------------*/
// module list
//
// Please keep the list of namespaces and the vector 'all' in sync.

namespace sp2 {
namespace python {
namespace fake_modules {

    static const char * PACKAGE = "_sp2_fake";

    // 1) Declarations of the fake modules, defined in separate files/namespaces
    namespace example { extern fake_module_t fake_module; }
    namespace mutation_helper { extern fake_module_t fake_module; }

    // 2) An iterable container over all fake modules.
    static const std::vector<fake_module_t*> all = {
        &example::fake_module,
        &mutation_helper::fake_module
    };

} // namespace fake_modules
} // namespace python
} // namespace sp2

#endif //SP2_PYTHON_FAKE_MODULES_HPP
