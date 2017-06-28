#ifndef SP2_BOND_CONTROL_T_HPP
#define SP2_BOND_CONTROL_T_HPP

/// \file fbc.hpp
/// \brief Main header for a fast bond control class scpp::fbc::bond_control_t

#include <vector>

#include "common/graph/ud_graph_t.hpp"
#include "common/neighbor/cell_array_t.hpp"

namespace sp2 {
/// fast bond control namespace
namespace fbc {

graph::ud_graph_t get_bond_graph(std::vector<double> positions,
    double lattice[3][3], double bond_max);

/// fast bond controller class
class bond_control_t
{
private:
    int na,                 ///< number of atoms from the last update() call
        nb;                 ///< number of bonds from the last update() call

    bool bond_lock;         ///< toggle whether or not bonds can be created/destroyed

    double cell_min,        ///< the minimum distance used to create cells (with an orthogonal basis, the cells are cubes with this side length)
        cell_padding;    ///< the portion of the cell_min value used to update only atoms that have moved more than (this amount / 2, cell_padding < cell_min)

    /// periodic cell array
    cell_array_t cell_array;

    std::vector<double> last_positions, ///< positions from the last update() call
        last_lattice_positions;         ///< lattice coordinate positions from the last update call

    /// full bond list (nb x nb)
    std::vector<std::vector<bool> > full_bond_list;
    std::vector<double> delta_list, ///< (nb x 3) matrix of bond components (x, y, z)
        length_list;                ///< (nb x 1) vector of bond lengths
    std::vector<int> bond_ids,      ///< (nb x 1) vector of bond ids
        sister_ids,                 ///< (nb x 1) vector for "sister" bond ids
        offsets;                    ///< ((na + 1) x 1) vector for bond id offsets for atoms

    int n_periodic;                     ///< number of (original) lattice vectors that had nonzero lengths, aka the number of periodic directions
    double orig_lattice[3][3];          ///< original lattice input

    double lattice[3][3],               ///< processed lattice vectors
        inv_lattice[3][3];           ///< inverse lattice matrix

    double rotation[3][3];              ///< rotation transformation calculated when processing lattice vectors

    double transformation[3][3],        ///< transformation matrix from normal to lattice coordinates
        inv_transformation[3][3];    ///< transformation matrix from lattice to normal coordinates

    /// bonds to be excluded from length and delta calculations
    std::vector<std::pair<int, int> > excluded_bonds;
public:
    /// constructor, class needs to be initialized with init() before use
    bond_control_t() : na(0), nb(0),
        bond_lock(false), cell_min(0), cell_padding(0) {}

    bond_control_t(const double input_lattice[3][3],
        double bond_max, double recalc_threshold) : bond_control_t()
    {
        init(input_lattice, bond_max, recalc_threshold);
    }

    /// \brief initialize the object
    /// \param input_lattice int[3][3] input lattice vectors (rows), leave 0 for non-periodic directions
    /// \param bond_max double bonds are guaranteed to exist if atoms are within this distance
    /// \param recalc_threshold distance that an atom must move before its bonds are recalculated, very large values will reduce performance
    void init(const double input_lattice[3][3],
        double bond_max, double recalc_threshold);

    /// \brief set the lattice vectors
    /// \param input_lattice int[3][3] input lattice vectors (rows), leave 0 for non-periodic directions
    void set_lattice(const double input_lattice[3][3]);

    /// \brief update function, calculates the bond list/information for the input
    /// \param positions vector<double>& input atom positions (na x 3)
    /// \exception atom_position_err if there is an issue with atom positions
    void update(const std::vector<double> &positions);

    /// \brief wrap input positions to within the current periodic cell
    /// \param positions vector<double>& input atom positions (na x 3)
    void wrap_positions(std::vector<double> &positions);

    /// get number of bonds
    int get_nb() const {return nb;}
    /// get number of atoms
    int get_na() const {return na;}
    /// get number of cells containing atoms
    int get_nc() const {return cell_array.n_active;}
    /// get number of periodic directions
    int get_periodicity() const {return n_periodic;}

    /// \brief calculate the components of the vector that connects the inputs, takes into account periodicity
    /// \param pos_a double* pointer to the coordinates of the first atom
    /// \param pos_b double* pointer to the coordinates of the second atom
    /// \param output double* output pointer for the vector components (pos_b - pos_a)
    void calc_delta(const double *pos_a, const double *pos_b, double *output) const;

    /// \brief calculate all delta vectors between atoms from the last update
    /// \param deltas vector<double>& output vector
    /// \param working vector<double>& working memory vector
    void calc_delta_all(std::vector<double> &deltas, std::vector<double> &working);

    /// \brief calculate all delta vectors between atoms as defined by the input graph and positions
    /// \param input_graph const util::graph_t& input undirected graph
    /// \param input_positions vector<double>& input atom positions
    /// \return vector<double> bond deltas the same order as in the input graph, in normal coordinates
    std::vector<double> calc_delta_graph(const graph::ud_graph_t &input_graph,
        const std::vector<double> &input_positions);

    /// \brief calculate the distance between the inputs, takes into account periodicity
    /// \param pos_a double* pointer to the coordinates of the first atom
    /// \param pos_b double* pointer to the coordinates of the second atom
    /// \return distance between pos_a and pos_b
    double calc_distance(const double *pos_a, const double *pos_b) const;

    /// transform input positions into lattice coordinates
    void transform_to_lattice(const std::vector<double> &positions,
        std::vector<double> &lattice_positions) const;
    ///  transform input lattice coordinates to regular coordinates
    void transform_from_lattice(const std::vector<double> &lattice_positions,
        std::vector<double> &positions) const;

    /// get bond delta vector
    const std::vector<double> &get_bond_deltas() const {return delta_list;}
    /// get bond length vector
    const std::vector<double> &get_bond_lengths() const {return length_list;}

    /// get bond id vector
    const std::vector<int> &get_bond_ids() const {return bond_ids;}
    /// get bond offset vector
    const std::vector<int> &get_bond_offsets() const {return offsets;}
    /// get sister bond id vector
    const std::vector<int> &get_sister_bonds() const {return sister_ids;}

    graph::ud_graph_t get_graph() const {
        return graph::ud_graph_t(offsets, bond_ids);}

    /// lock bonds and do not allow them to be created or destroyed when update() is called
    void lock_bonds() {bond_lock = true;}
    /// unlock bonds and allow them to be created or destroyed when update() is called
    void unlock_bonds() {bond_lock = false;}
    /// directly set the object's bonds
    void set_bonds(const std::vector<int> &input_offsets,
        const std::vector<int> &input_bond_ids);

    /// exclude the creation of a bond between atoms id_a and id_b
    void set_excluded_bond(int id_a, int id_b);
    /// allow creation of a bond between atoms id_a and id_b
    void remove_excluded_bond(int id_a, int id_b);
    /// clear all excluded bonds
    void clear_excluded_bonds();

private:
    /// disallow copy constructor
    bond_control_t(const bond_control_t &);

    /// special function for small systems where "boxing" may mess up
    void update_small();

    /// function to resize vectors and values dependent on the number of atoms in the system
    void resize_members(int num_atoms);

    /// check changes in position
    void mark_movement(const std::vector<double> &positions,
        std::vector<double> &delta);
    /// get bond list and delta values from input positions
    void get_bonds(std::vector<double> &lattice_positions,
        std::vector<double> &temp_deltas,
        std::vector<std::pair<int, int> > &temp_bonds);
    /// process bond list based on delta values
    void process_bonds(std::vector<double> &temp_deltas,
        std::vector<std::pair<int, int> > &temp_bonds);
};

} // namespace fbc
} // namespace sp2

#endif // SP2_BOND_CONTROL_T_HPP
