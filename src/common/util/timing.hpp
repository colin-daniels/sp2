//
// Created by cc on 1/28/17.
//

#ifndef SP2_CL_TIMER_T_HPP
#define SP2_CL_TIMER_T_HPP

#include <time.h>
#include <functional>

namespace sp2 {

double benchmark(unsigned int iters, std::function<void(void)> func);

/// a small timer class.
class cl_timer_t
{
private:
    timespec m_start, ///< start time
        m_end;   ///< end time
public:
    cl_timer_t() : m_start(), m_end()
    {
        m_start.tv_nsec = 0;
        m_end.tv_nsec = 0;
    };

/// start the timer and clear the last run
    void start()
    {
        clock_gettime(CLOCK_REALTIME, &m_start);
        m_end.tv_sec = 0;
    };

/// stop the timer and save the elapsed time
    void stop()
    { clock_gettime(CLOCK_REALTIME, &m_end); };

/// get seconds since timer start or last duration if the timer is stopped
    double duration()
    {
        if (m_end.tv_sec == 0)
        {
            timespec temp;
            clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &temp);
            return (temp.tv_sec - m_start.tv_sec) +
                   1e-9 * (temp.tv_nsec - m_start.tv_nsec);
        }
        else
            return (m_end.tv_sec - m_start.tv_sec) +
                   1e-9 * (m_end.tv_nsec - m_start.tv_nsec);
    };
};

} // namespace sp2

#endif // SP2_CL_TIMER_T_HPP
