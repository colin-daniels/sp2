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
#include "common/minimize/metropolis.hpp" // re-export

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
[[nodiscard]] std::vector<double> acgsd(
    scalar_fn_t objective_fn, vector_fn_t gradient_fn,
    std::vector<double> initial_position,
    const acgsd_settings_t &settings = {});
/// \brief Non-linear conjugate gradient minimization routine, ACGSD \cite andrei2006c .
[[nodiscard]] std::vector<double> acgsd(diff_fn_t objective_fn,
    std::vector<double> initial_position,
    const acgsd_settings_t &settings = {});

class uphill_linesearch: public std::exception
{
    virtual const char* what() const noexcept
    { return "Uphill linesearch"; }
};

/// \brief Cubic/quadratic interpolation linesearch.
/// \param objective_fn : Input function to be minimized.
/// \param slope_fn : The objective function's associated scalar gradient function.
/// \param alpha : The initial stepsize for the linesearch.
/// \returns The best value of alpha found before exit.
double linesearch(const ls_settings_t &settings,
    double alpha, oned_fn_t objective_fn, oned_fn_t slope_fn);

/// \brief Cubic/quadratic interpolation linesearch.
double linesearch(const ls_settings_t &settings,
    double alpha, diff1d_fn_t diff_fn);

/// \brief Linear conjugate gradient for sparse matrices (solves Ax = b).
/// \param matrix_fn : Function that returns the result of Ax for some x
/// \param b : The vector b in Ax = b.
/// \param tolerance : Exit condition, ||residual||_{maximum norm} < tolerance to exit.
/// \returns The vector x within tolerance.
[[nodiscard]] std::vector<double> linear_cg(vector_fn_t matrix_fn,
    std::vector<double> b, double tolerance = 1e-12);

[[nodiscard]] std::vector<double> fire(diff_fn_t grad_fn,
    double mass, std::vector<double> initial_position,
    const fire_settings_t &settings = {});

// PSO currently cannot be run without MPI
#ifdef SP2_ENABLE_MPI

/// \brief Adaptive particle swarm optimization \cite zhan2009adaptive .
/// \param settings : Input settings.
/// \param objective_fn : The input objective function.
/// \param comm : The mpi communicator on which the pso will run.
/// \returns The best solution found before exit.
[[nodiscard]] std::vector<double> adaptive_pso(scalar_fn_t objective_fn,
    const pso_settings_t &settings, boost::mpi::communicator comm);

/// function should evaluate current position and set value for particle
class particle_t;
using pso_update_fn_t = std::function<void(particle_t&)>;

/// specialization for update functions that need to know more or modify
/// the calling particle in APSO
[[nodiscard]] std::vector<double> adaptive_pso(pso_update_fn_t update_fn,
    const pso_settings_t &settings, boost::mpi::communicator comm);

#endif // SP2_ENABLE_MPI

} // namespace minimize
} // namespace sp2

#endif // SP2_MINIMIZE_HPP
