#ifndef SP2_METROPOLIS_ENUMS_HPP
#define SP2_METROPOLIS_ENUMS_HPP

// This belongs to metropolis.  However, apparently it is illegal to specialize
// sp2::enum_map from within sp2::minimize, so it has its own file to avoid the
// whole namespace song and dance.

#include "common/enums.hpp"
#include "common/math/vec3_util.hpp"
#include "common/python/types/as_ndarray.hpp"

namespace sp2 {

/// types of changes that can be made to a structure
enum class structural_mutation_type : int
{
    INVALID = 0,
    CART_COORDS = 1,
    FRAC_COORDS = 2,
    LATTICE = 3
};

template<>
constexpr enum_map_t<structural_mutation_type>
    enum_map<structural_mutation_type> = {
    {structural_mutation_type::CART_COORDS, "carts"},
    {structural_mutation_type::FRAC_COORDS, "fracs"},
    {structural_mutation_type::LATTICE,     "lattice"}
};

/// type of objective function to use for a minimization method that
///  doesn't explicitly use force
enum class objective_type : int
{
    NONE      = 0,
    POTENTIAL = 1,
    GRAD_NORM = 2,
    GRAD_MAX  = 3
};

template<>
constexpr enum_map_t<objective_type> enum_map<objective_type> = {
    {objective_type::POTENTIAL, "potential"},
    {objective_type::GRAD_NORM, "grad_norm"},
    {objective_type::GRAD_MAX,  "grad_max"}
};

/// result of a mutation function on a structure.
struct structural_mutation_t
{
    structural_mutation_type type;
    as_ndarray_t<double> data;
};

static double compute_objective(objective_type type,
    const std::pair<double, std::vector<double>> &diff)
{
    switch (type)
    {
    case objective_type::POTENTIAL:
        return diff.first;

    case objective_type::GRAD_NORM:
    {
        double sum = 0;
        for (auto x : diff.second)
            sum += x * x;
        return sqrt(sum);
    }

    case objective_type::GRAD_MAX:
    {
        double best = 0;
        for (auto &x : dtov3(diff.second))
            best = std::max(best, x.mag());
        return best;
    }

    default:
        throw std::logic_error("unrecognized objective type");
    }
};

} // namespace sp2

#endif // SP2_METROPOLIS_ENUMS_HPP
