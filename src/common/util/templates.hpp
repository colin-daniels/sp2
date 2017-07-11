#ifndef SP2_TEMPLATES_HPP
#define SP2_TEMPLATES_HPP

#include <vector>
#include <algorithm>
#include <type_traits>
#include <functional>

namespace sp2 {

template<class T>
struct ni;

template<class F>
class scope_guard_t
{
private:
    bool exec = false;
    F destructor;

public:
    scope_guard_t(F destructor)
        : exec{true}, destructor{destructor} {}

    scope_guard_t(scope_guard_t &&other)
        : exec{true}, destructor{std::move(other.destructor)}
    {
        other.exec = false;
    }

    ~scope_guard_t()
    {
        if (exec)
            destructor();
    }
};

template<class F>
inline auto scope_guard(F &&func) { return scope_guard_t<F>{func}; }

template<class IterBegin, class IterEnd>
struct range_wrapper_t
{
    IterBegin it_begin;
    IterEnd it_end;

    decltype(auto) begin() const { return it_begin; }
    decltype(auto) begin()       { return it_begin; }

    decltype(auto) end() const { return it_end; }
    decltype(auto) end()       { return it_end; }
};

template<class IterBegin, class IterEnd>
inline auto make_range(IterBegin &&begin, IterEnd &&end)
{
    return range_wrapper_t<IterBegin, IterEnd>{
        std::forward<IterBegin>(begin),
        std::forward<IterEnd>(end)
    };
}

template<class T = std::size_t>
struct id_iter_t
{
    T current,
        max;

    constexpr T operator*() const { return current; }
    constexpr id_iter_t& operator++() { ++current; return *this; }
    constexpr bool operator!=(const id_iter_t&) const {
        return current != max; }
};

template<class T>
auto id_range(const T &container)
{
    using size_type = decltype(container.size());
    return make_range(
        id_iter_t<size_type>{0, container.size()},
        id_iter_t<size_type>{}
    );
}

inline auto id_range(std::size_t begin, std::size_t end)
{
    return make_range(
        id_iter_t<>{begin, end},
        id_iter_t<>{}
    );
}

template<class T, std::size_t N>
auto id_range(const T (&)[N])
{
    return make_range(
        id_iter_t<std::size_t>{0, N},
        id_iter_t<std::size_t>{}
    );
}

// TODO: fix reverse_range (what's wrong with it?)
template<class T>
auto reverse_range(const T &container)
{
    return make_range(
        container.rbegin(),
        container.rend()
    );
}

/// makes a vector of length n via std::generate and the provided generator
template<class T, class Generator>
std::vector<T> make_vector(size_t n, Generator gen,
    std::result_of_t<Generator()>* = nullptr)
{
    auto result = std::vector<T>();
    std::generate_n(std::back_inserter(result), n, gen);
    return result;
}

template<class Generator>
auto make_vector(size_t n, Generator gen,
    std::result_of_t<Generator()>* = nullptr)
{
    return make_vector<std::decay_t<
        std::result_of_t<Generator()>>>(n, gen);
}

/// makes a vector of length n via std::generate, the provided generator,
/// and additionally passes the current vector index to the generator
/// for each element generated
template<class T, class Generator>
std::vector<T> make_vector(size_t n, Generator gen,
    std::result_of_t<Generator(std::size_t)>* = nullptr)
{
    std::size_t idx = 0;
    return make_vector<T>(n, [&](){
        return gen(idx++);
    });
}

template<class Generator>
auto make_vector(size_t n, Generator gen,
    std::result_of_t<Generator(std::size_t)>* = nullptr)
{
    return make_vector<std::decay_t<
        std::result_of_t<Generator(std::size_t)>>>(n, gen);
}

/// return the current index of the element that would be in
/// the nth position if the input range was sorted
template<class RandomIt>
std::size_t nth_index(RandomIt first, RandomIt last, std::size_t n)
{
    auto len = std::distance(first, last);

    std::vector<std::size_t> indices(len);
    for (std::size_t i = 0; i < indices.size(); ++i)
        indices[i] = i;

    auto nth = indices.begin() + n;
    std::nth_element(indices.begin(), nth, indices.end(),
        [&](size_t a, size_t b) {
            return *(first + a) < *(first + b);
        });

    return *nth;
}

} // namespace sp2

#endif // SP2_TEMPLATES_HPP
