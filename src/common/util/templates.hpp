#ifndef SP2_TEMPLATES_HPP
#define SP2_TEMPLATES_HPP

#include <vector>
#include <algorithm>
#include <type_traits>
#include <functional>

namespace sp2 {

template<class T>
struct not_implemented;

template<class T>
struct range_wrap_t
{
    T const it_begin;
    T const it_end;

    range_wrap_t(T input_begin, T input_end) :
        it_begin(input_begin), it_end(input_end) {}

    const T& begin() const {return it_begin;}
    const T& end() const {return it_end;}
};

template<class T>
range_wrap_t<T> range(T const first, T const last) {
    return {first, last};}

struct id_iter_t
{
    std::size_t i;

    id_iter_t(std::size_t inp) : i(inp) {}

    std::size_t operator*() const {return i;}
    id_iter_t& operator++() {++i; return *this;}
    bool operator!=(const id_iter_t &other) const {
        return i != other.i;}
};

template<class T>
range_wrap_t<id_iter_t> id_range(const T &container)
{
    return range_wrap_t<id_iter_t>(
        id_iter_t(0),
        id_iter_t(container.size())
    );
}

template<class T, std::size_t N>
range_wrap_t<id_iter_t> id_range(const T (&)[N])
{
    return range_wrap_t<id_iter_t>(
        id_iter_t(0),
        id_iter_t(N)
    );
}

// TODO: fix reverse_range
template<class T>
auto reverse_range(const T &container)
{
    return range_wrap_t<decltype(container.rbegin())>{
        container.rbegin(),
        container.rend()
    };
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
