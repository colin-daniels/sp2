//
// Created by cc on 9/7/16.
//

#ifndef SP2_STRUCTURE_T_HPP
#define SP2_STRUCTURE_T_HPP

#include "common/json/json_serializable_t.hpp"
#include "common/enums.hpp"

#include <vector>
#include <boost/mpi/communicator.hpp>
#include <common/math/vec3_t.hpp>
#include <common/math/mat3x3_t.hpp>

namespace sp2 {

/// structural information, atom positions/types/lattice vectors
struct structure_t : public io::json_serializable_t
{
    /// lattice vectors (each row is a vector), a zero lattice vector
    /// is taken to mean that the direction is non-periodic
    double lattice[3][3] = {};

    /// space group name (if applicable), international notation
    std::string space_group;
    int n_symm;

    std::vector<atom_type> types;   ///< N x 1 vector of atom types (see: enums.hpp)
    std::vector<vec3_t> positions;  ///< N x 1 vector of cartesian positions

    structure_t() = default;
    structure_t(const double lattice_in[3][3],
        const std::vector<atom_type> &types_in,
        const std::vector<double> &positions_in);

    structure_t(const double lattice_in[3][3],
        const std::vector<atom_type> &types_in,
        const std::vector<vec3_t> &positions_in);

    /// Rotate the structure, leaving the lattice fixed.
    void rotate(vec3_t axis, double theta);
    /// Apply a general cartesian transformation matrix, modifying the positions
    /// (but leaving the lattice fixed).
    void transform(const mat3x3_t &transformation);

    bool serialize(Json::Value &output) const;
    bool deserialize(const Json::Value &input);

    Json::Value serialize_lattice() const;
    bool deserialize_lattice(const Json::Value &input);

    /// MPI broadcast member function.
    void bcast(const boost::mpi::communicator &comm, int root);

    /* ------------------------------------------------------------ */
    // Working with fractional positions

    /// Compute fractional positions for the atoms.
    ///
    /// Caveats:
    /// - Non-periodic axes are not currently supported.
    /// - The positions may not necessarily lay within [0.0, 1.0).
    ///   (See reduced_fractional_positions()).
    std::vector<vec3_t> fractional_positions() const;

    /// Compute fractional positions for the atoms, reduced into
    /// the primitive cell.
    ///
    /// Caveats:
    /// - Non-periodic axes are not currently supported.
    std::vector<vec3_t> reduced_fractional_positions() const;

    /// Set fractional positions for the atoms.
    ///
    /// Caveats:
    /// - Non-periodic axes are not currently supported.
    void set_fractional_positions(const std::vector<vec3_t> &positions);

    /// Reduce positions into the confines of the lattice primitive cell.
    ///
    /// Be aware that this may affect scaling operations, as atoms are moved to
    /// equivalent locations closer to the origin; this may or may not be
    /// desired, depending on the application!
    ///
    /// Caveats:
    /// - Non-periodic axes are not currently supported.
    void reduce_positions();

    /* ------------------------------------------------------------ */
    // Retrieving lattice properties

    /// Compute the unit cell volume as a positive quantity.
    ///
    /// Caveats:
    /// - Output is unspecified in the presence of non-periodic axes.
    constexpr double get_volume() const
    {
        return abs(mat3x3_t(lattice).determinant());
    }

    /// Compute the length of each lattice vector.
    ///
    /// Caveats:
    /// - Output is unspecified in the presence of non-periodic axes.
    vec3_t get_lengths() const;

    /* ------------------------------------------------------------ */
    // Rescaling (modifying the lattice at constant fractional positions)

    /// Scale the structure by a relative factor along each lattice vector.
    ///
    /// Caveats:
    /// - Non-periodic axes are not currently supported.
    /// - Consider if reduce_positions() should be performed first.
    void scale(vec3_t factors);

    /// Take a linear combination of the lattice vectors, preserving
    /// fractional coordinates.  This is just the straightforward
    /// generalization of scale(vec3_t) beyond diagonal scale matrices.
    ///
    /// Caveats:
    /// - Non-periodic axes are not supported.
    /// - Consider if reduce_positions() should be performed first.
    void scale(mat3x3_t factor);

    /// Individually scale each lattice vector to a specified length.
    ///
    /// Caveats:
    /// - Non-periodic axes are not currently supported.
    /// - Consider if reduce_positions() should be performed first.
    void rescale(vec3_t new_lengths);

    /// Set a new lattice, preserving fractional coordinates.
    ///
    /// Caveats:
    /// - Non-periodic axes are not currently supported.
    /// - Consider if reduce_positions() should be performed first.
    void rescale(mat3x3_t new_lattice);

    // These have been decided against since they're ambiguous; are we
    // scaling lengths or volumes? (only volume makes sense for rescale(double),
    // but it still seems like scale(double) could be easily misused)

    // void scale(double factor);
    // void rescale(double volume);
};

void sort_structure_types(structure_t &structure);

} // namespace sp2

#endif // SP2_STRUCTURE_T_HPP
