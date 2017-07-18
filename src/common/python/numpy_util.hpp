#ifndef SP2_NUMPY_UTIL_HPP
#define SP2_NUMPY_UTIL_HPP

#include <vector>

/// Represents a numpy ndarray.
///
/// This mostly just exists as a target for Python (de)serialization.
template <typename T>
class ndarray_serialize_t {
    // The dimensions of the array.  Slowest index first (C order).
    std::vector<std::size_t> _shape;
    // The elements of the array, in C order.
    std::vector<T> _data;

    // Invariant:  _data.size() = product(_shape)

public:

    // Make a 1D array of length 0.
    ndarray_serialize_t() : _shape{0} {}

    // Make a 0D array (scalar).
    ndarray_serialize_t(T value) : _data{value} {}

    // Make a 1D array.
    ndarray_serialize_t(std::vector<T> v) : _data(v), _shape{v.size()} {}

    // Make an arbitrary ndarray.
    ndarray_serialize_t(std::vector<T> data, std::vector<std::size_t> shape) {
        auto product = 1;
        for (auto x : shape)
            product *= x;

        if (product != data.size()) {
            throw std::length_error("Length of input data inconsistent with shape");
        }

        _data = move(data);
        _shape = move(shape);
    }

    const std::vector<T>& data() const { return _data; }

    size_t size() const { return _data.size(); }
    size_t ndim() const { return _shape.size(); }
    const std::vector<size_t>& shape() const { return _shape; }
};


#endif //SP2_NUMPY_UTIL_HPP
