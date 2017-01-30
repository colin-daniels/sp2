#ifndef SP2_PERIODIC_CELL_T_HPP
#define SP2_PERIODIC_CELL_T_HPP

/// \file periodic_cell_t.hpp
/// \brief Definition header for the scpp::fbc::periodic_cell_t class

#include <vector>
#include <string>
#include <stdexcept>

namespace sp2 {
namespace fbc {

/// single periodic cell class
class periodic_cell_t
{
private:
    int id,                         ///< cell id
        time_empty;                 ///< number of updates that this cell has had zero atoms
    periodic_cell_t *neighbors[13], ///< forward neighbor links
        *p_neighbors[13];           ///< backward neighbor links
    std::vector<int> atom_ids,      ///< ids of atoms within the cell
        to_check;                   ///< atoms that are to be checked for new bonds
    std::vector<int> comp;          ///< atoms that the to_check atoms will be checked against

    friend class cell_array_t;
    friend class bond_control_t;

    /// maximum number of atoms within a single unit cell
    static const unsigned int atom_max = 256;
public:
    /// constructor, input is the cell id
    periodic_cell_t(int i) noexcept;
    /// destructor, automatically unlinks the cell
    ~periodic_cell_t() {unlink_all();}

    /// link one cell to another with a specific neighbor position
    void link(periodic_cell_t* other, int i) {
        neighbors[i] = other; other->p_neighbors[i] = this;}
    /// unlink from all cells
    void unlink_all();

    /// \brief add an atom to the cell
    /// \exception cell_num_atoms_err if the number of atoms in a
    /// single cell has exceeded the maximum allowed
    void add_atom(const int id, std::vector<int> &last_icx);
    /// remove an atom from the cell
    void remove_atom(const int id, std::vector<int> &last_icx);
    /// get the number of atoms in the cell
    int  size() {return atom_ids.size();}

    /// clear to_check and increment time_empty if applicable
    void update_clear();
    /// add an atom to the to_check vector
    void add_check(const int i) {to_check.push_back(i);}
    /// update the get_comp vector
    void get_comp();
};

inline periodic_cell_t::periodic_cell_t(int i) noexcept :
    id(i), time_empty(0)
{
    for (int j = 0; j < 13; ++j)
    {
        neighbors[j] = nullptr;
        p_neighbors[j] = nullptr;
    }
}

inline void periodic_cell_t::unlink_all()
{
    for (int i = 0; i < 13; ++i)
    {
        if (neighbors[i] != nullptr)
            neighbors[i]->p_neighbors[i] = nullptr;
        if (p_neighbors[i] != nullptr)
            p_neighbors[i]->neighbors[i] = nullptr;

        neighbors[i] = nullptr;
        p_neighbors[i] = nullptr;
    }
}

inline void periodic_cell_t::add_atom(const int index,
    std::vector<int> &last_icx)
{
    if (atom_ids.size() > atom_max)
        throw std::runtime_error(
            "number of atoms in periodic_cell_t  exceeded the maximum of "
            + std::to_string(atom_max));

    last_icx[index] = atom_ids.size();
    atom_ids.push_back(index);
}

inline void periodic_cell_t::remove_atom(const int index,
    std::vector<int> &last_icx)
{
    last_icx[atom_ids.back()] = last_icx[index];
    atom_ids[last_icx[index]] = atom_ids.back();
    atom_ids.pop_back();
}

inline void periodic_cell_t::update_clear()
{
    to_check.clear();
    if (atom_ids.empty())
        ++time_empty;
    else
        time_empty = 0;
}

inline void periodic_cell_t::get_comp()
{
    comp.clear();
    for (int i = 0; i < 13; ++i)
        if (neighbors[i] != nullptr)
            for (unsigned int j = 0; j < neighbors[i]->to_check.size(); ++j)
                comp.push_back(neighbors[i]->to_check[j]);
}

} // namespace fbc
} // namespace sp2

#endif // SP2_PERIODIC_CELL_T_HPP
