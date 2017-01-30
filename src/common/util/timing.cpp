#include "common/util/timing.hpp"

double sp2::benchmark(unsigned int iters, std::function<void(void)> func)
{
    sp2::cl_timer_t ct;

    ct.start();
    for (unsigned int i = 0; i < iters; ++i)
        func();
    ct.stop();

    return ct.duration();
}