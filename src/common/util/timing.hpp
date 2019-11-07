#ifndef SP2_CL_TIMER_T_HPP
#define SP2_CL_TIMER_T_HPP

#include <functional>
#include <chrono>

namespace sp2 {

double benchmark(unsigned int iters, std::function<void(void)> func);

/// a small timer class.
class cl_timer_t
{
public:
    /// start the timer and clear the last run
    void start()
    {
        time_start = clock_type::now();
        running = true;
    }

    /// stop the timer and save the elapsed time (returned by duration)
    void stop()
    {
        time_stop = clock_type::now();
        running = false;
    }

    /// get seconds since timer start or last duration if the timer is stopped
    double duration()
    {
        const auto end = (running ? clock_type::now() : time_stop);
        return std::chrono::duration_cast<std::chrono::duration<double>>(
            end - time_start).count();
    }

    bool is_running() const { return running; }

private:
    using clock_type = std::conditional_t<
        // use high resolution clock if it's steady
        std::chrono::high_resolution_clock::is_steady,
        std::chrono::high_resolution_clock,
        // otherwise fall back to steady clock
        std::chrono::steady_clock
    >;

    bool running = false;
    std::chrono::time_point<clock_type> time_start,
        time_stop;
};

} // namespace sp2

#endif // SP2_CL_TIMER_T_HPP
