#include <Python.h> // Must be first include

#include "bindings.hpp"

#include "modules/fake_modules.hpp"
#include "modules/ext_modules.hpp"
#include "common/python/types/py_object_t.hpp"

//-----------------------------

using namespace std;
using namespace sp2::python;
namespace sp2 {
namespace python {

py_object_t sp2::python::run_phonopy::make_extra_kw(
    std::vector<size_t> sc_to_prim)
{
    py_ref_t py_sc_map = [&] {
        auto list = py_from(sc_to_prim);
        auto &module = fake_modules::mutation_helper::fake_module.module;
        auto klass = getattr(module, "supercell_index_mapper");

        auto args = scope(Py_BuildValue("(O)", list.inner().raw()));
        throw_on_py_err("Exception constructing python args tuple.");


        auto po = py_object_t(
            std::shared_ptr<py_ref_t>(&klass, [](auto) {})
        );

#warning FIXME
        //return call_callable(klass, just_args(args.dup()));
        return list.inner();
#warning FIXME
    }();

    py_ref_t kw = [&] {
        auto kw = scope(Py_BuildValue("{sO}",
            "supercell", py_sc_map.raw()
        ));
        throw_on_py_err("Exception constructing python kw dict.");
        return kw.move();
    }();

    return opaque(kw);
}

} // namespace sp2
} // namespace python
