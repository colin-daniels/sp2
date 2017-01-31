//
// Created by cc on 9/7/16.
//

#ifndef SP2_RANDOM_HPP
#define SP2_RANDOM_HPP

#include <random>

#ifdef SP2_ENABLE_MPI
#include <boost/mpi/communicator.hpp>
#endif // SP2_ENABLE_MPI

namespace sp2 {
namespace util {

/// Simple random number generator convenience class.
class rng_t
{
private:
    /// Underlying random number generator.
    std::mt19937_64 gen;
    // Distributions as member vars to avoid reconstructing every
    // time rand() is called.
    std::uniform_int_distribution<int> int_dist;
    std::uniform_real_distribution<double> real_dist;

public:
    /// Seed the generator with a given seed sequence.
    void seed(std::seed_seq&& seq);
    /// Seed the generator with a single integer.
    void seed(uint64_t num);
    /// Seed the generator using random_device (/dev/urandom probably).
    void seed_random(std::size_t len = 2);

    /// Random double [0, 1).
    double rand();
    /// Random double [low, high).
    double rand(double low, double high);
    /// Random integer [low, high).
    int rand(int low, int high);
    /// Get underlying pseudorandom number generator object.
    auto& get_gen() {return gen;}

#ifdef SP2_ENABLE_MPI
    /// MPI broadcast member function.
    void bcast(const boost::mpi::communicator &comm, int root);
#endif // SP2_ENABLE_MPI

};

inline double rng_t::rand()
{
    return std::generate_canonical<double, 64>(gen);
}

inline double rng_t::rand(double low, double high)
{
    using param_t = decltype(real_dist)::param_type;
    return real_dist(gen, param_t{low, high});
}

inline int rng_t::rand(int low, int high)
{
    using param_t = decltype(int_dist)::param_type;
    return int_dist(gen, param_t{low, high - 1});
}

} // namespace util
} // namespace sp2


#endif //SP2_RANDOM_HPP
