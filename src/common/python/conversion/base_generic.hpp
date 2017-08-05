#ifndef SP2_PYTHON_CONVERSION_BASE_GENERIC_HPP
#define SP2_PYTHON_CONVERSION_BASE_GENERIC_HPP

#include "base_generic_fwd.hpp"
#include "base_monomorphic_fwd.hpp"

// base_generic.hpp:
//
//  Provides '{to,from}_python' on generic types for 'py_object_t'.
//  These implementations are forced to rely on the small interface provided
//  by 'py_object_t', and are probably laughably inefficient compared to
//  the same operations on 'py_ref_t'.

//      !!!!!!!!!!!!!!!!!!!
//        !!  CAUTION  !!
//
// 'base_generic.hpp' and 'base_generic_raw.hpp' both independently implement
// the same functions on different types.  Never edit one file without editing
// the other.
//
// Unfortunately there is currently no viable abstraction that would enable
// a single implementation to be written for both types.
//
//        !!  CAUTION  !!
//      !!!!!!!!!!!!!!!!!!!

#include "concrete.hpp"
#include "numpy.hpp"

#include <type_traits>

namespace sp2 {
namespace python {

/* -------------------------------------------------------------------------- */
// Generic conversions of:
//
//  - std::vector ---> list
//  - std::vector <--- sequence or iterable

template<typename T>
bool to_python(const std::vector<T> &vec, py_object_t &py_list)
{
    if (py_list)
        throw std::logic_error(
            "attempted to serialize into occupied py_object_t");

    py_list = py_import("builtins").getattr("list").call();
    auto append = py_list.getattr("append");

    for (const auto &x : vec)
    {
        py_object_t py_item;
        if (!to_python(x, py_item))
        {
            py_list.destroy();
            return false;
        }

        append.call({py_item});
    }

    return true;
}

template<typename T>
bool from_python(const py_object_t &o, std::vector<T> &vec)
{
    vec.clear();

    // We want to accept any *iterable*, but not *iterators*;
    // Consuming one-time use iterators would be uncool.
    if (o.hasattr("__next__"))
    {
        std::cerr << "error: refusing to deserialize a one-time use iterator. "
                  << "Please use e.g. list() to strictly evaluate it."
                  << std::endl;
        return false;
    }

    // iterable into sequence
    auto py_sequence = py_import("builtins").getattr("list").call({o});

    // use sequence protocol
    std::size_t size;
    if (!from_python(py_sequence.getattr("__len__").call(), size))
        return false;

    vec.reserve(size);
    auto getitem = py_sequence.getattr("__getitem__");
    for (std::size_t i = 0; i < size; i++)
    {
        py_object_t py_i;
        if (!to_python(i, py_i))
        {
            vec.clear();
            return false;
        }

        auto py_item = getitem.call({py_i});

        T cxx_item;
        if (!from_python(py_item, cxx_item))
        {
            vec.clear();
            return false;
        }

        vec.push_back(std::move(cxx_item));
    }
    return true;
}

/* -------------------------------------------------------------------------- */
// Generic conversions of:
//
//  - std::tuple ---> list (FIXME: justify why not tuple)
//  - std::tuple <--- sequence or iterable

template<std::size_t Idx, class ...Ts>
constexpr bool _from_python_rec(Ts &&...)
{ return true; }

template<std::size_t Idx = 0, class ...Ts,
    class = std::enable_if_t<(Idx < sizeof...(Ts))>>
bool _from_python_rec(std::vector<py_object_t> &pys, std::tuple<Ts...> &cs)
{
    return from_python(pys[Idx], std::get<Idx>(cs)) &&
           _from_python_rec<Idx + 1>(pys, cs);
}

template<typename... Ts>
bool from_python(const py_object_t &py, std::tuple<Ts...> &cs)
{
    std::vector<py_object_t> pys;
    if (!from_python(py, pys))
        return false;

    if (pys.size() != sizeof...(Ts))
    {
        std::cerr << "expected python sequence of length " << sizeof...(Ts)
                  << ", got length " << pys.size() << std::endl;
        return false;
    }

    return _from_python_rec(pys, cs);
}

// -----

template<std::size_t Idx>
constexpr bool _to_python_rec(...)
{ return true; }

template<std::size_t Idx = 0, class ...Ts,
    class = std::enable_if_t<(Idx < sizeof...(Ts))>>
bool _to_python_rec(const std::tuple<Ts...> &cs, std::vector<py_object_t> &pys)
{
    return to_python(std::get<Idx>(cs), pys[Idx]) &&
           _to_python_rec<Idx + 1>(cs, pys);
}

template<typename... Ts>
bool to_python(const std::tuple<Ts...> &cs, py_object_t &py)
{
    std::vector<py_object_t> pys(sizeof...(Ts));

    if (!_to_python_rec(cs, pys))
    {
        // destruct any python references that were created prior to the error
        pys.resize(0);
    }

    return to_python(pys, py);
}

/* -------------------------------------------------------------------------- */
} // namespace python
} // namespace sp2

#endif // SP2_PYTHON_CONVERSION_BASE_GENERIC_HPP
