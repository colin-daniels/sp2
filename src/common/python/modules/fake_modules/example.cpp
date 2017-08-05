#include "Python.h" // must be first import.

#include "common/python/modules/fake_modules.hpp"

namespace sp2 {
namespace python {
namespace fake_modules {
namespace example {

fake_module_t fake_module({
    // NOTE: Pretend this says "example"
    //       This is mangled only because there is also an ext_module example
    "_example_fake",
    "\n\
def test():                                                                  \n\
    print('you did it')                                                      \n\
"});

} // namespace example
} // namespace fake_modules
} // namespace python
} // namespace sp2
