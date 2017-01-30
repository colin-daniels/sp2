#ifndef SP2_CELL_ARRAY_T_HPP
#define SP2_CELL_ARRAY_T_HPP

/// \file cell_array_t.hpp
/// \brief Definition header for the scpp::fbc::cell_array_t class

// periodic_cell_t definition
#include "common/neighbor/periodic_cell_t.hpp"

#include <vector>

namespace sp2 {
namespace fbc {

/// periodic cell array class
class cell_array_t
{
private:
    int n_cells,    ///< number of periodic cells in the array
        n_active,   ///< number of periodic cells in the array
        n_periodic; ///< number of periodic directions

    /// cell 'lifetime', how many updates until a cell with no atoms in it is automatically destroyed
    int cell_lifetime;

    int low_bound[3],   ///< minimum lattice coordinates
        high_bound[3],  ///< maximum lattice coordinates
        range[3];       ///< number of cells in each lattice direction

    /// pointer vector containing entries for all cells within the array's limits defined by low_bound and high_bound
    std::vector<periodic_cell_t*> cells_full;

    std::vector<int> last_indices,  ///< atom->cell indices from the last time update_cells() was called
        last_icx;                   ///< atom->cell->atom_ids indices from the last time update_cells() was called

    friend class bond_control_t;
public:
    /// normal constructor, class must still be initialized before use
    cell_array_t();
    /// destructor, deletes all periodic_cell_t objects
    ~cell_array_t();

    /// \brief initialize
    /// \param input_low int[3] lower index bound (inclusive) for each lattice vector direction
    /// \param input_high int[3] upper index bound (exclusive) for each lattice vector direction
    /// \param num_periodic int how many lattice directions are periodic (value from 0 to 3)
    void init(int input_low[3], int input_high[3], int num_periodic);

    /// \brief update cells
    /// \param lattice_positions vector<double>& array of the lattice coordinates of all the atoms
    /// \param delta vector<double>& vector designating (via value != 0) atoms that will be added to each cell's bond update list
    /// \exception cell_index_err if there is an issue with calculated indices
    /// \exception cell_num_atoms_err if the number of atoms in a single cell has exceeded the maximum allowed
    void update_cells(const std::vector<double> &lattice_positions,
        std::vector<double> &delta);

    /// \brief resize operation
    /// \param input_low int[3] minimum lattice coordinates
    /// \param input_high int[3] maximum lattice coordinates
    void resize_array(int input_low[3], int input_high[3]);

    /// get number of cells (including empty ones)
    int size() const {return cells_full.size();}

    /// get cell ranges
    void get_range(int r[3]) const {
        for (int i = 0; i < 3; ++i) {r[i] = range[i];};}

private:
    /// disallow copy constructor
    cell_array_t(const cell_array_t &);

    /// calculate cell indices for the input positions
    void get_indices(const std::vector<double> &lattice_positions,
        std::vector<int> &indices);

    /// construct periodic cell at the location
    void add_cell(int x, int y, int z);
    /// destroy the periodic cell at the location
    void remove_cell(int x, int y, int z);

    /// link all cells to their neighbors
    void link_all();
    /// unlink all cells from their neighbors
    void unlink_all();

    /// link a single cell to all of its neighbors
    void link_cell(int id);
    /// link a single cell to all of its neighbors
    void link_cell(int x, int y, int z);

    /// convert cell xyz indices to a cell index
    int index(int x, int y, int z);
};

} // namespace fbc
} // namespace sp2

#endif // SP2_CELL_ARRAY_T_HPP
