#ifndef SP2_NC_UTILITY_FUNCTIONS_HPP
#define SP2_NC_UTILITY_FUNCTIONS_HPP

/// \file utility_functions.hpp
/// \brief Utility functions for scpp::fbc

namespace sp2 {
namespace fbc {

/// \brief sanitize the lattice input
/// \param lattice double[3][3] input lattice vectors in matrix form
/// \param mod_lattice double[3][3] the modified lattice (rotated, and populated if less than 3 inputs were given)
/// \param transformation double[3][3] the 3x3 rotation matrix that was applied to the lattice input
/// \return int number of input lattice vectors that were nonzero (number of periodic directions)
/// \exception lattice_err if the processed lattice has a zero determinant
int process_lattice(const double lattice[3][3],
    double mod_lattice[3][3], double transformation[3][3]);

void invert_3x3(const double input[3][3], double inverse[3][3]);

void array_rot(const int axis, double *input,
    double *output, const double theta);

} // namespace fbc
} // namespace sp2

#endif // SP2_NC_UTILITY_FUNCTIONS_HPP
