#include "common/vec3_t.hpp"
#include "common/util/random.hpp"

using namespace std;
using namespace sp2;

template<class R>
void randomize_vec(double *pos, R rng_fn)
{
    // x: between -1 and 1
    pos[0] = 1.0 - 2.0 * rng_fn();

    // y: cos(theta) * sqrt(1 - x^2)
    // z: sin(theta) * sqrt(1 - x^2)
    const double two_pi = 2 * atan2(0, -1);
    sincos(
            two_pi * rng_fn(),
            &pos[1], // sin
            &pos[2]  // cos
    );

    pos[1] *= sqrt(1 - pos[0] * pos[0]);
    pos[2] *= sqrt(1 - pos[0] * pos[0]);
}

vec3_t& vec3_t::randomize()
{
    randomize_vec(pos, []() {
        return rand() / (RAND_MAX + 1.0);
    });

    return *this;
}

vec3_t& vec3_t::randomize(util::rng_t &rng)
{
    randomize_vec(pos, [&]() {
        return rng.rand(0.0, 1.0);
    });

    return *this;
}

vec3_t sp2::unit_normal(const vec3_t &v)
{
    vec3_t normal_vec(0, 0, 0);

    int n_non_zero = 0;
    for (auto d : v)
        n_non_zero += (d != 0);

    if (n_non_zero == 0)
        normal_vec[0] = 1;
    else if (n_non_zero == 1 || n_non_zero == 2)
    {
        for (int i = 0; i < 3; ++i)
        {
            if (v[i] == 0)
            {
                normal_vec[i] = 1;
                break;
            }
        }
    }
    else
    {
        normal_vec[0] =  v[1];
        normal_vec[1] = -v[0];
        normal_vec.normalize();
    }

    return normal_vec;
}

std::vector<double> sp2::v3tod(const std::vector<vec3_t> &input)
{
    vector<double> output;
    output.reserve(input.size() * 3);
    for (auto &v : input)
        for (auto &d : v)
            output.push_back(d);

    return output;
}

std::vector<sp2::vec3_t> sp2::dtov3(const std::vector<double> &input)
{
    vector<vec3_t> output;
    output.reserve(input.size() / 3);
    for (size_t i = 0; i < input.size(); i += 3)
        output.emplace_back(&input[i]);

    return output;
}
