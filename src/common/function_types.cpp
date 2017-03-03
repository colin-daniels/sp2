#include "common/function_types.hpp"
#include "common/math/blas.hpp"

using namespace std;
using namespace sp2;

oned_fn_t sp2::to_one_dim(scalar_fn_t input_fn,
                           const std::vector<double> &initial_pos,
                           const std::vector<double> &direction)
{
    return [=](double x) {
        auto pos = initial_pos;
        vaxpy(x, direction, pos);
        return input_fn(pos);
    };
}

oned_fn_t sp2::to_one_dim(scalar_fn_t input_fn,
                           const std::vector<double> &direction)
{
    return [=](double x) {
        auto pos = direction;
        vscal(x, pos);
        return input_fn(pos);
    };
}

diff1d_fn_t sp2::to_one_dim(diff_fn_t input_fn,
                             const std::vector<double> &initial_pos,
                             const std::vector<double> &direction)
{
    return [=](double x) {
        auto pos = initial_pos;
        vaxpy(x, direction, pos);

        auto result = input_fn(pos);
        return make_pair(result.first, vdot(result.second, direction));
    };
}

diff1d_fn_t sp2::to_one_dim(diff_fn_t input_fn,
                             const std::vector<double> &direction)
{
    return [=](double x) {
        auto pos = direction;
        vscal(x, pos);

        auto result = input_fn(pos);
        return make_pair(result.first, vdot(result.second, direction));
    };
}
