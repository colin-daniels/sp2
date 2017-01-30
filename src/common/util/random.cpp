//
// Created by cc on 9/7/16.
//

#include "common/util/random.hpp"

#ifdef SP2_USE_MPI
#include <boost/mpi.hpp>

namespace mpi = boost::mpi;
#endif // SP2_USE_MPI

using namespace std;
using namespace sp2;

void util::rng_t::seed(std::seed_seq&& seq)
{
    gen.seed(seq);
}

void util::rng_t::seed(uint64_t num)
{
    gen.seed(num);
}

void util::rng_t::seed_random(std::size_t len)
{
    std::random_device dev;
    std::vector<decltype(dev())> seeds;
    for (auto i = 0U; i < len; ++i)
        seeds.push_back(dev());

    seed(std::seed_seq(seeds.begin(), seeds.end()));
}

#ifdef SP2_USE_MPI
void util::rng_t::bcast(const mpi::communicator &comm, int root)
{
    // serialize to a string
    stringstream oss;
    oss << gen;

    // broadcast the state
    auto state = oss.str();
    mpi::broadcast(comm, state, root);

    // deserialize from the newly-broadcasted state string
    stringstream iss(state);
    iss >> gen;
}
#endif // SP2_USE_MPI