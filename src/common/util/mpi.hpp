#ifndef SP2_MPI_HPP
#define SP2_MPI_HPP

#include <utility>
#include <type_traits>
#include <boost/mpi/collectives.hpp>

////////////////////////////////////////////////////////////////////////////////
// Miscellaneous MPI utility functions                                        //
////////////////////////////////////////////////////////////////////////////////

namespace boost {
namespace mpi {

/// Minimum location collective operation.
template <class T> using min_loc = minimum<std::pair<T, int>>;
/// Maximum loaction collective operation.
template <class T> using max_loc = maximum<std::pair<T, int>>;

// Overloads for min/max loc operations.
template<typename T>
int all_reduce(const communicator &comm, const T &in_value, min_loc<T> op);
template<typename T>
int all_reduce(const communicator &comm, const T &in_value, max_loc<T> op);

} // namespace mpi
} // namespace boost


////////////////////////////////////////////////////////////////////////////////
// Implementations                                                            //
////////////////////////////////////////////////////////////////////////////////

namespace boost {
namespace mpi {

/// Convenience broadcast template for user types with a
/// bcast(const communicator&, int) member function.
template<typename T>
void broadcast(const communicator &comm, T &value, int root,
    std::enable_if_t<std::is_function<typename T::bcast>::value>* = nullptr)
{
    value.bcast(comm, root);
}

template<typename T, typename U>
void broadcast(const communicator &comm, std::pair<T, U> &value, int root)
{
    broadcast(comm, value.first,  root);
    broadcast(comm, value.second, root);
}

template<typename T>
int all_reduce(const communicator &comm, const T &in_value, min_loc<T> op)
{
    std::pair<T, int> out_pair;
    all_reduce(comm, std::make_pair(in_value, comm.rank()), out_pair, op);
    return out_pair.second;
}

template<typename T>
int all_reduce(const communicator &comm, const T &in_value, max_loc<T> op)
{
    std::pair<T, int> out_pair;
    all_reduce(comm, std::make_pair(in_value, comm.rank()), out_pair, op);
    return out_pair.second;
}

} // namespace mpi
} // namespace boost

#endif // SP2_MPI_HPP
