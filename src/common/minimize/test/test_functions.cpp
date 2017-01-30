#include "test_functions.hpp"

using namespace std;
using namespace sp2;

namespace {
const size_t dim = 25;
}

const minimize::test_fn_t minimize::sphere_fn{
    // dimensions
    dim,
    // minimum value
    0.0,
    // minimum location
    vector<double>(dim, 0.0),
    // lower bound on search domain
    vector<double>(dim, -5),
    // upper bound on search domain
    vector<double>(dim, +5),
    // evaluation function
    [](const vector<double> &pos) {
        double value = 0;
        vector<double> gradient;
        for (size_t i = 0; i < dim; ++i)
        {
            value += pos[i] * pos[i];
            gradient.push_back(2 * pos[i]);
        }
        return make_pair(value, gradient);
    }
};

const minimize::test_fn_t minimize::sum_squares_fn{
    // dimensions
    dim,
    // minimum value
    0.0,
    // minimum location
    vector<double>(dim, 0.0),
    // lower bound on search domain
    vector<double>(dim, -10),
    // upper bound on search domain
    vector<double>(dim, +10),
    // evaluation function
    [](const vector<double> &pos) {
        double value = 0;
        vector<double> gradient;
        for (size_t i = 0; i < dim; ++i)
        {
            value += (i + 1) * pos[i] * pos[i];
            gradient.push_back(2 * (i + 1) * pos[i]);
        }
        return make_pair(value, gradient);
    }
};

const minimize::test_fn_t minimize::styblinski_tang_fn{
    // dimensions
    dim,
    // minimum value
    -39.16616570377142 * dim,
    // minimum location
    vector<double>(dim, -2.903534018185960),
    // lower bound on search domain
    vector<double>(dim, -5),
    // upper bound on search domain
    vector<double>(dim, +5),
    // evaluation function
    [](const vector<double> &pos) {
        double value = 0;
        vector<double> gradient;
        for (auto x : pos)
        {
            value += 0.5 * (x * x * x * x - 16 * x * x + 5 * x);
            gradient.push_back(2 * x * x * x - 16 * x + 2.5);
        }
        return make_pair(value, gradient);
    }
};

const minimize::test_fn_t minimize::rosenbrock_fn{
    // dimensions
    dim,
    // minimum value
    0.0,
    // minimum location
    vector<double>(dim, 1),
    // lower bound on search domain
    vector<double>(dim, -5),
    // upper bound on search domain
    vector<double>(dim, +10),
    // evaluation function
    [](const vector<double> &pos) {
        double value = 0;
        vector<double> gradient(dim, 0);
        for (size_t i = 0; i + 1 < dim; ++i)
        {
            double v1 = (pos[i + 1] - pos[i] * pos[i]),
                v2 = (pos[i] - 1);
            value += 100 * v1 * v1 + v2 * v2;
            gradient[i] += -400 * pos[i] * v1 + 2 * v2;
            gradient[i + 1] += 200 * v1;
        }

        return make_pair(value, gradient);
    }
};

const minimize::test_fn_t minimize::perm_fn{
    // dimensions
    dim,
    // minimum value
    0.0,
    // minimum location
    [](){
        vector<double> min;
        for (size_t i = 0; i < dim; ++i)
            min.push_back(1.0 / (i + 1));
        return min;
    }(),
    // lower bound on search domain
    vector<double>(25, -25),
    // upper bound on search domain
    vector<double>(25, +25),
    // evaluation function
    [](const vector<double> &pos) {
        double value = 0;
        vector<double> gradient(dim, 0);

        for (size_t i = 0; i < dim; ++i)
        {
            double sum = 0;
            for (size_t j = 0; j < dim; ++j)
                sum += (j + 1.0) * (pow(pos[j], i + 1)
                                    - 1.0 / pow(j + 1.0, i + 1));

            value += sum * sum;
            for (size_t j = 0; j < dim; ++j)
                gradient[j] += 2 * (i + 1) * pow(pos[j], i) * (j + 1.0) * sum;
        }

        return make_pair(value, gradient);
    }
};

