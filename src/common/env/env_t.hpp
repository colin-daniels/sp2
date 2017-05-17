#ifndef SP2_ENV_T_HPP
#define SP2_ENV_T_HPP

#include <randutils.hpp>
#include <mutex>
#include <boost/mpi/communicator.hpp>
#include <sstream>
#include <fstream>
#include <common/util/mpi.hpp>
#include "run/program_args_t.hpp"

template<typename T>
using mixer_rep = std::array<typename T::result_type, T::size()>;

auto get_randutils_param(randutils::seed_seq_fe256 &sseq)
{
    mixer_rep<randutils::seed_seq_fe256> state;
    sseq.param(state.begin());

    return state;
}

std::string to_string(randutils::seed_seq_fe256 &sseq)
{
    auto state = get_randutils_param(sseq);

    std::ostringstream oss;
    oss << '{';
    for (std::size_t i = 0; i < state.size(); ++i)
    {
        oss << "\n\t0x"
            << std::hex << std::uppercase
            << state[i] << 'u';
        if (i + 1 != state.size())
            oss << ',';
    }
    oss << "\n}\n";

    return oss.str();
}

namespace sp2 {

class logger_t
{
public:

    template<typename ...Args>
    void log_info(Args &&...args)
    {
        log_impl(std::cout, std::forward<Args>(args)...);
    }

private:
//    std::ofstream stream_dev,
//        stream_info,
//        stream_warning,
//        stream_error,
//        stream_fatal;

    void log_impl(std::ostream &) {}

    template<typename T, typename ...Args>
    void log_impl(std::ostream &stream, T&& object, Args &&...args)
    {
        stream << std::forward<T>(object);
        log_impl(stream, std::forward<Args>(args)...);
    }
};

class env_t
{
public:
    env_t() = delete;
    env_t(const env_t&) = delete;
    env_t& operator=(const env_t&) = delete;

    env_t(env_t&&) = default;
    env_t& operator=(env_t&&) = default;

//    env_t(int & argc, char **& argv);

private:
    boost::mpi::communicator comm;
    logger_t logger;
    randutils::mt19937_rng rng;


    env_t(const boost::mpi::communicator &comm_in, const logger_t &logger_in,
        randutils::mt19937_rng &&rng_in);

    friend std::vector<env_t> construct_thread_environments(unsigned int,
        const boost::mpi::communicator &, logger_t,
        randutils::seed_seq_fe256 &&);
};


//    randutils::mt19937_rng get_new_rng()
//    {
//        // construct and return a new rng for the requesting thread
//        randutils::mt19937_rng rng(seed_seq);
//        // note we need to stir seed_seq since otherwise we will return the same
//        // constructed rng every time this function is called
//        seed_seq.stir();
//        return rng;
//    }

std::vector<env_t> construct_thread_environments(
    unsigned int n_threads,
    const boost::mpi::communicator &comm,
    logger_t logger,
    randutils::seed_seq_fe256 &&seeds
)
{
    if (n_threads == 0)
        throw std::invalid_argument("Number of threads cannot be zero.");

    std::vector<env_t> environments;
    environments.reserve(n_threads);

    logger.log_info("initial seed: ", to_string(seeds));

    for (unsigned int i = 0; i < n_threads; ++i)
    {
        environments.emplace_back(
            env_t(comm, logger, randutils::mt19937_rng(seeds))
        );

        seeds.stir();
    }

    return environments;
}


//std::vector<env_t> construct_thread_environments(
//    unsigned int n_threads = std::thread::hardware_concurrency(),
//    const boost::mpi::communicator &comm = {},
//    logger_t logger = {})
//{
//    randutils::seed_seq_fe256 seeds = {0};
//    if (comm.rank() == 0)
//        seeds = std::move(randutils::auto_seed_256{}.base());
//
//    auto param = get_randutils_param(seeds);
//    boost::mpi::broadcast(comm, param, 0);
//
//    std::vector<mixer_rep<decltype(seeds)>> seed_param;
//    if (comm.rank() == 0)
//    {
//        randutils::auto_seed_256 auto_seeds;
//        seed_param.push_back(get_randutils_param(auto_seeds.base()));
//
//        for (int i = 1; i < comm.size(); ++i)
//        {
//
//        }
//
//    }
////    get_randutils_param
//
//    return construct_thread_environments(
//        n_threads, comm, logger,
//        std::move(seeds));
//}
//
//std::vector<randutils::mt19937_rng> get_seeded_rngs(
//    const boost::mpi::communicator &comm,
//    int n_rngs
//)
//{
//    return {};
//}
//
//std::vector<randutils::mt19937_rng> get_seeded_rngs(
//    const boost::mpi::communicator &comm,
//    int n_rngs,
//    randutils::seed_seq_fe256 &seed
//)
//{
//    using seed_vec_t = std::vector<mixer_rep<decltype(seed)>>;
//
//    std::vector<seed_vec_t> seeds_to_scatter;
//    if (comm.rank() == 0)
//    {
//
//        std::vector<> seed_param;
//    }
//
//    return {};
//}


} // namespace sp2


#endif // SP2_ENV_T_HPP
