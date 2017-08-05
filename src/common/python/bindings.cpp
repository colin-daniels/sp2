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

py_object_t structural_metropolis::make_param_pack(vector<double> carts,
    const double lattice[3][3], vector<double> force)
{
    // interpret as 3N cartesian coords
    size_t width = 3;
    if (carts.size() % width != 0)
        throw logic_error("vector not divisible by width");

    size_t height = carts.size() / width;

    auto py_carts = py_from(as_ndarray(carts, {height, width}));

    auto c_lattice = vector<double>(&lattice[0][0], &lattice[0][0] + 9);
    auto py_lattice = py_from(as_ndarray(c_lattice, {3, 3}));

    auto py_force = py_from(as_ndarray(force, {height, width}));

    py_ref_t kw = [&] {
        auto kw = scope(Py_BuildValue("{sOsOsO}",
            "carts", py_carts.inner().raw(),
            "lattice", py_lattice.inner().raw(),
            "force", py_force.inner().raw()
        ));
        throw_on_py_err("Exception constructing python kw dict.");
        return move(kw);
    }();

    return opaque(kw);
}

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
