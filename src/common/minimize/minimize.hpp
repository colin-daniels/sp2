#ifndef SP2_MINIMIZE_HPP
#define SP2_MINIMIZE_HPP

/// \file minimize.hpp
/// \brief Main header for minimization functions.

#include <vector>
#include <functional>

#ifdef SP2_ENABLE_MPI
#include <boost/mpi.hpp>
#endif // SP2_ENABLE_MPI

#include "common/minimize/settings.hpp"
#include "common/function_types.hpp"

namespace sp2 {
/// Minimization method namespace.
namespace minimize {

////////////////////////////////////////////////////////////////////////////////
// Minimization Methods                                                       //
////////////////////////////////////////////////////////////////////////////////

/// \brief Non-linear conjugate gradient minimization routine, ACGSD \cite andrei2006c .
/// \param objective_fn : Input function to be minimized.
/// \param gradient_fn : The objective function's associated vector gradient function.
/// \param initial_position : The initial position to start at.
/// \returns The best solution found before exit.
std::vector<double> acgsd(scalar_fn_t objective_fn, vector_fn_t gradient_fn,
    std::vector<double> initial_position,
    const acgsd_settings_t &settings = {});
/// \brief Non-linear conjugate gradient minimization routine, ACGSD \cite andrei2006c .
std::vector<double> acgsd(diff_fn_t objective_fn,
    std::vector<double> initial_position,
    const acgsd_settings_t &settings = {});

/// \brief Cubic/quadratic interpolation linesearch.
/// \param objective_fn : Input function to be minimized.
/// \param slope_fn : The objective function's associated scalar gradient function.
/// \param alpha : The initial stepsize for the linesearch.
/// \returns The best value of alpha found before exit.
double linesearch(oned_fn_t objective_fn, oned_fn_t slope_fn, double alpha);
/// \brief Cubic/quadratic interpolation linesearch.
double linesearch(diff1d_fn_t objective_fn, double alpha);

/// \brief Linear conjugate gradient for sparse matrices (solves Ax = b).
/// \param matrix_fn : Function that returns the result of Ax for some x
/// \param b : The vector b in Ax = b.
/// \param tolerance : Exit condition, ||residual||_{maximum norm} < tolerance to exit.
/// \returns The vector x within tolerance.
std::vector<double> linear_cg(vector_fn_t matrix_fn,
    std::vector<double> b, double tolerance = 1e-12);

std::vector<double> fire(diff_fn_t grad_fn,
    double mass, std::vector<double> initial_position,
    const fire_settings_t &settings = {});

/// \brief Minimization algorithm using random mutations
/// \param objective_fn : Function to be minimized
/// \param mutation_fn : Function that applies a random mutation in-place
/// \returns The position of the minimal value found
std::vector<double> metropolis(
    scalar_fn_t objective_fn,
    mutation_fn_t mutation_fn,
    std::vector<double> initial_position,
    const metropolis_settings_t &settings = {});

// PSO currently cannot be run without MPI
#ifdef SP2_ENABLE_MPI

/// \brief Adaptive particle swarm optimization \cite zhan2009adaptive .
/// \param settings : Input settings.
/// \param objective_fn : The input objective function.
/// \param comm : The mpi communicator on which the pso will run.
/// \returns The best solution found before exit.
std::vector<double> adaptive_pso(scalar_fn_t objective_fn,
    const pso_settings_t &settings, boost::mpi::communicator comm);

/// function should evaluate current position and set value for particle
class particle_t;
using pso_update_fn_t = std::function<void(particle_t&)>;

/// specialization for update functions that need to know more or modify
/// the calling particle in APSO
std::vector<double> adaptive_pso(pso_update_fn_t update_fn,
    const pso_settings_t &settings, boost::mpi::communicator comm);

#endif // SP2_ENABLE_MPI

} // namespace minimize
} // namespace sp2

#endif // SP2_MINIMIZE_HPP
