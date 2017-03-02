#include "common/vec3_t.hpp"
#include "common/util/random.hpp"

using namespace std;
using namespace sp2;

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
