//
// Created by cc on 1/29/17.
//

#ifndef SP2_UTIL_HPP
#define SP2_UTIL_HPP

#include <type_traits>
#include <cinttypes>
#include <cstring>

#include <templ/map.hpp>

/// copy const specifier from T to U
template<class T, class U>
using copy_const_t = std::conditional_t<
    std::is_const<T>::value, std::add_const<U>, U>;

/// copy volatile specifier from T to U
template<class T, class U>
using copy_volatile_t = std::conditional_t<
    std::is_volatile<T>::value, std::add_volatile<U>, U>;

/// copy cv specifiers from type T to U
template<class T, class U>
using copy_cv_t = copy_const_t<T, copy_volatile_t<T, U>>;

/// get unsigned char pointer to examine the input's object representation
template<class T>
auto get_uchar_ptr(T& input)
{
    using ptr_type = copy_cv_t<T, unsigned char*>;
    return reinterpret_cast<ptr_type>(&input);
}



/// closest (larger or equal) power of 2 (should only be used in constexpr
/// expressions due to being non-optimal)
constexpr std::size_t ceil_pow2(std::size_t input)
{
    std::size_t power = 1;
    while (power < input)
        power <<= 1;

    return power;
}



/// convert obeject to unsigned integer, preserving its bitwise representation
template<
    class T,
    // select int type based on size of input type
    class UIntType = templ::switch_t<
        // we round up to the closest power of 2 make sure T fits in the output
        templ::ic_size_t<ceil_pow2(sizeof(T))>,
        // type map from size to unsigned integer type
        templ::keyval<templ::ic_size_t< sizeof(uint8_t)>,  uint8_t>,
        templ::keyval<templ::ic_size_t<sizeof(uint16_t)>, uint16_t>,
        templ::keyval<templ::ic_size_t<sizeof(uint32_t)>, uint32_t>,
        templ::keyval<templ::ic_size_t<sizeof(uint64_t)>, uint64_t>,
        // default case
        templ::keyval<templ::ic_size_t<sizeof(T)>, uint64_t>
    >,
    // if the default type was selected we need to make sure the size of our
    // result type is large enough
    class = std::enable_if_t<(sizeof(T) <= sizeof(UIntType))>
>
inline UIntType as_uint_t(T input)
{
    // additional assert just in case
    static_assert(sizeof(T) <= sizeof(UIntType), "");

    UIntType result{};
    std::memcpy(&result, &input, sizeof(input));

    return result;
}

/// convert object from unsigned integer, preserving bitwise representation
template<
    class T, class UIntType,
    class = std::enable_if_t<(sizeof(T) <= sizeof(UIntType))>
>
inline T from_uint_t(UIntType input)
{
    // additional assert just in case
    static_assert(sizeof(T) <= sizeof(UIntType), "");

    T result{};
    std::memcpy(&result, &input, sizeof(result));
    return result;
}


#endif // SP2_UTIL_HPP
