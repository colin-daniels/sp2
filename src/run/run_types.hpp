#ifndef SP2_RUN_TYPES_HPP
#define SP2_RUN_TYPES_HPP

#include <string>
#include <json/json.h>
#include <mpi.h>

#include "common/structure_t.hpp"
#include "common/enums.hpp"
#include "run/run_settings_t.hpp"

namespace sp2 {
////////////////////////////////////////////////////////////////////////////////
// Main run types                                                             //
////////////////////////////////////////////////////////////////////////////////

/// execute energy minimization run/structural relaxation
int run_minimize(const run_settings_t &config, MPI_Comm comm);

/// execute a run of ATAC
int run_atac(const run_settings_t &config, MPI_Comm comm);

/// run the symmetry structure search algorithm
int run_symm(const run_settings_t &config, MPI_Comm comm);

/// run phonopy to calculate raman spectra
int run_phonopy(const run_settings_t &config, MPI_Comm comm);

////////////////////////////////////////////////////////////////////////////////
// Utility + common functions                                                 //
////////////////////////////////////////////////////////////////////////////////

/// output configuration defaults to the specified file
int generate_defaults(std::string filename);

} // namespace sp2

#endif // SP2_RUN_TYPES_HPP
