//
// Created by cc on 9/7/16.
//

#ifndef SP2_STRUCTURE_T_HPP
#define SP2_STRUCTURE_T_HPP

#include "common/json/json_serializable_t.hpp"
#include "common/enums.hpp"

#include <vector>
#include <boost/mpi/communicator.hpp>

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
    std::vector<double> positions;  ///< N x 3 matrix of atom positions,
    /// e.g. (x1, y1, z1, x2, y2, z2, ...)

    structure_t() = default;
    structure_t(const double lattice_in[3][3],
        const std::vector<atom_type> &types_in,
        const std::vector<double> &positions_in);

    bool serialize(Json::Value &output) const;
    bool deserialize(const Json::Value &input);

    Json::Value serialize_lattice() const;
    bool deserialize_lattice(const Json::Value &input);

    /// MPI broadcast member function.
    void bcast(const boost::mpi::communicator &comm, int root);
};

void sort_structure_types(structure_t &structure);

} // namespace sp2

#endif // SP2_STRUCTURE_T_HPP
