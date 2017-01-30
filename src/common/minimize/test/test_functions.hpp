#ifndef SP2_TEST_FUNCTIONS_HPP
#define SP2_TEST_FUNCTIONS_HPP

#include <vector>
#include "common/minimize/minimize.hpp"

namespace sp2 {
namespace minimize {

struct test_fn_t
{
    size_t dimensions;
    double min_value;
    std::vector<double> minimum,
        min_bound,
        max_bound;
    diff_fn_t function;

    auto operator()(const std::vector<double> &pos) const {
        return function(pos);}
};

// optimization test functions
extern const test_fn_t sphere_fn;
extern const test_fn_t sum_squares_fn;
extern const test_fn_t styblinski_tang_fn;
extern const test_fn_t rosenbrock_fn;
extern const test_fn_t perm_fn;


} // namespace minimize
} // namespace sp2


#endif //SP2_TEST_FUNCTIONS_HPP
