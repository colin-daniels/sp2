#ifndef SP2_TEST_HPP_HPP
#define SP2_TEST_HPP_HPP

#include <utility>
#include <type_traits>
#include <cmath>
#include <tuple>
#include <array>

namespace ptnl {

////////////////////////////////////////////////////////////////////////////////
// AIREBO potential calculation helper functions                              //
////////////////////////////////////////////////////////////////////////////////

enum class atom_types : int
{
    HYDROGEN = 1,
    CARBON = 2
};

enum class bond_types : int
{
    HH = 0,
    HC = 1,
    CC = 3
};

/// get atom-atom bond type from constituent atom types a and b
constexpr bond_types get_bond_type(atom_types a, atom_types b);


/// switching function, t < 0 = 1 and t > 1 = 0
constexpr auto switching_fn(const double t);

/// switching function, t = 0 -> 1, t = 1 -> 0, assumes t is between [0:1]
constexpr auto switching_fn_nocheck(const double t);


template<bond_types bond_type>
constexpr auto cutoff_function(const double r);

template<bond_types bond_type>
constexpr auto attractive_potential(const double r);

template<>
constexpr auto attractive_potential<bond_types::CC>(const double r);

template<bond_types bond_type>
constexpr auto repulsive_potential(const double r);


template<bond_types bond_type>
constexpr auto lennard_jones_potential(const double r);

template<bond_types bond_type>
constexpr auto torsion_potential(const double cos_omega);


/// carbon G(theta)
constexpr auto gtheta_c(const double cos_theta, const double N_t);

/// hydrogen G(theta)
constexpr auto gtheta_h(const double cos_theta);

////////////////////////////////////////////////////////////////////////////////
// helper lookup-table exponential function                                   //
////////////////////////////////////////////////////////////////////////////////

namespace detail {

#include "generated/exp_table.hpp"

/// exponential lookup table function, valid in the range [-9.9:0]
constexpr double lk_exp(double x)
{
    constexpr auto step = exp_table_step,
        inv_step = 1.0 / step;

    const double dst = x - static_cast<int>(inv_step * x) * step - step / 2;
    return (1.0 + dst + (dst * dst) / 2)
           * exp_table[static_cast<int>(inv_step * x)];
}

} // namespace detail

////////////////////////////////////////////////////////////////////////////////
// utility aliases                                                            //
////////////////////////////////////////////////////////////////////////////////

template<class T>
constexpr decltype(auto) val(T &&arg)
{
    return std::get<0>(std::forward<T>(arg));
}

template<class T>
constexpr decltype(auto) deriv(T &&arg)
{
    return std::get<1>(std::forward<T>(arg));
}

template<std::size_t N, class T>
constexpr decltype(auto) deriv_n(T &&arg)
{
    return std::get<N>(std::forward<T>(arg));
}

//constexpr auto value = 0U;
//constexpr auto deriv = 1U;

////////////////////////////////////////////////////////////////////////////////
// atom/bond types and related functions                                      //
////////////////////////////////////////////////////////////////////////////////

constexpr bond_types get_bond_type(atom_types a, atom_types b)
{
    // instead of using if/else we just do some arithmetic
    return static_cast<bond_types>(
        static_cast<int>(a) * static_cast<int>(b) - 1
    );
}

static_assert(get_bond_type(
    atom_types::CARBON, atom_types::CARBON) == bond_types::CC, "");
static_assert(get_bond_type(
    atom_types::HYDROGEN, atom_types::CARBON) == bond_types::HC, "");
static_assert(get_bond_type(
    atom_types::HYDROGEN, atom_types::HYDROGEN) == bond_types::HH, "");

namespace detail {

template<bond_types type, typename U>
constexpr U get_bond_constant(bond_types input_type, U value)
{
    return input_type == type ?
           value : U{};
}

template<bond_types type, typename U, typename ...Args>
constexpr U get_bond_constant(bond_types input_type, U value, Args &&...args)
{
    return input_type == type ?
           value : get_bond_constant<type>(std::forward<Args>(args)...);
}

} // namespace detail

////////////////////////////////////////////////////////////////////////////////
// potential calculation functions                                            //
////////////////////////////////////////////////////////////////////////////////

constexpr auto switching_fn_nocheck(const double t)
{
    // the switching function is the Hermite basis function h_00
    return std::make_tuple(
        // S(t)
        (1 + 2 * t) * (t - 1) * (t - 1),
        // d/dt S(t)
        6 * t * (t - 1)
    );
}

/// switching function, t < 0 = 1 and t > 1 = 0
constexpr auto switching_fn(const double t)
{
    if (t < 0)
        return std::make_tuple(1.0, 0.0);
    else if (t > 1)
        return std::make_tuple(0.0, 0.0);
    else
        return switching_fn_nocheck(t);
}

template<bond_types bond_type>
constexpr auto lennard_jones_potential(const double r)
{
    constexpr double epsilon = detail::get_bond_constant<bond_type>(
        bond_types::HH, 0.00150,
        bond_types::HC, 0.00206, // sqrt(epsilon_CC * epsilon_HH)
        bond_types::CC, 0.00284
    );

    constexpr double sigma = detail::get_bond_constant<bond_type>(
        bond_types::HH, 2.650,
        bond_types::HC, 3.025, // 1/2 (sigma_CC + sigma_HH)
        bond_types::CC, 3.400
    );

////////////////////////////////////////////////////////////////////////////////
// actual calculation                                                         //
////////////////////////////////////////////////////////////////////////////////

    // (sigma / r)^3
    const double sr3 = (sigma / r) * (sigma / r) * (sigma / r);

    // V^LJ(r) = 4 epsilon [(sigma / r)^12 - (sigma / r)^6]
    // d/dr V^LJ(r) = -24 epsilon [2 (sigma / r)^12 - (sigma / r)^6] / r
    return std::make_tuple(
        // V^LJ(r)
        4 * epsilon * (sr3 * sr3 * sr3 - sr3 * sr3),
        // d/dr V^LJ(r)
        -12 * epsilon * (2 * sr3 * sr3 * sr3 - sr3 * sr3) / r
    );
}

template<bond_types bond_type>
constexpr auto torsion_potential(const double cos_omega)
{
    // epsilon depending on atoms involved in the bond
    constexpr double epsilon = detail::get_bond_constant<bond_type>(
        bond_types::HH, 0.1250, // HC-CH
        bond_types::HC, 0.1787, // CC-CH and HC-CC
        bond_types::CC, 0.3079  // CC-CC
    );

    // V^tors(omega) = epsilon [(256 / 405) cos^10(omega/2) - 1/10]
    //  or
    // V^tors(omega) = epsilon [(8 / 405) (cos(omega) + 1)^5 - 1/10]

    // 1/2 [cos(x) + 1] = cos^2(x/2)
    const double x = cos_omega + 1;

    return std::make_tuple(
        // V^tors(cos(omega))
        epsilon * ((8.0 / 405.0) * (x * x) * (x * x) * x - 0.1),
        // d/dcos(omega) V^tors(cos(omega))
        epsilon * (8.0 /  81.0) * (x * x) * (x * x)
    );
}

template<bond_types bond_type>
constexpr auto cutoff_function(const double r)
{
////////////////////////////////////////////////////////////////////////////////
// constants used for calculation, change depending on bond type              //
////////////////////////////////////////////////////////////////////////////////

    // minimum cutoff distance for bond types in Angstroms
    constexpr double r_min = detail::get_bond_constant<bond_type>(
        bond_types::HH, 1.1,
        bond_types::HC, 1.3,
        bond_types::CC, 1.7
    );

    // maximum cutoff distance for bond types in Angstroms
    constexpr double r_max = detail::get_bond_constant<bond_type>(
        bond_types::HH, 1.7,
        bond_types::HC, 1.8,
        bond_types::CC, 2.0
    );

    constexpr double r_range = r_max - r_min;

////////////////////////////////////////////////////////////////////////////////
// actual calculation                                                         //
////////////////////////////////////////////////////////////////////////////////

    if (r < r_min)
        return std::make_tuple(1.0, 0.0);
    else if (r > r_max)
        return std::make_tuple(0.0, 0.0);
    else
    {
        // the cutoff function f^c(r) is just the switching function
        // with t = (r - r_min) / (r_max - r_min)
        //
        // note we use the nocheck version because we already made sure
        // that the bond length is in range
        auto s = switching_fn_nocheck((r - r_min) / r_range);

        return std::make_tuple(
            // f^c(r) aka just the value of the switching function
            val(s),
            // d/dr f^c(r) aka ds/dt * dt/dr
            deriv(s) / r_range
        );
    }
}

template<bond_types bond_type>
constexpr auto attractive_potential(const double r)
{
////////////////////////////////////////////////////////////////////////////////
// constants used for calculation, change depending on bond type (CH and HH)  //
////////////////////////////////////////////////////////////////////////////////

    constexpr double B = detail::get_bond_constant<bond_type>(
        bond_types::HH, 29.6325930000,
        bond_types::HC, 32.3551866587
    );

    constexpr double beta = detail::get_bond_constant<bond_type>(
        bond_types::HH, 1.71589217000,
        bond_types::HC, 1.43445805925
    );

////////////////////////////////////////////////////////////////////////////////
// actual potential and potential derivative calculation                      //
////////////////////////////////////////////////////////////////////////////////

    // V^A(r) = B exp(-beta r)
    return std::make_tuple(
        // V^A(r)
        B * detail::lk_exp(-beta * r),
        // d/dr V^A(r)
        -beta * B * detail::lk_exp(-beta * r)
    );
}

template<>
constexpr auto attractive_potential<bond_types::CC>(const double r)
{
////////////////////////////////////////////////////////////////////////////////
// potential calculation, specialized for CC bonds                            //
////////////////////////////////////////////////////////////////////////////////

    // scaling constants in eV
    constexpr double B[3] = {
        12388.79197798,
        17.56740646509,
        30.71493208065
    };

    // exponential constants in Angstroms^-1
    constexpr double beta[3] = {
        4.72045231270,
        1.43321324990,
        1.38269125060
    };

    // for CC bonds we have a three term sum instead of one term like HH/HC
    // V^A(r) = \sum_{n=1..3} B_n exp(-beta_n r)
    return std::make_tuple(
        // V^A(r)
          B[0] * detail::lk_exp(-beta[0] * r)
        + B[1] * detail::lk_exp(-beta[1] * r)
        + B[2] * detail::lk_exp(-beta[2] * r),
        // d/dr V^A(r)
        - beta[0] * B[0] * detail::lk_exp(-beta[0] * r)
        - beta[1] * B[1] * detail::lk_exp(-beta[1] * r)
        - beta[2] * B[2] * detail::lk_exp(-beta[2] * r)
    );
}

template<bond_types bond_type>
constexpr auto repulsive_potential(const double r)
{
////////////////////////////////////////////////////////////////////////////////
// constants used for calculation, change depending on bond type              //
////////////////////////////////////////////////////////////////////////////////

    // exponential constants by bond type in Angstroms^-1
    constexpr double alpha = detail::get_bond_constant<bond_type>(
        bond_types::HH, 3.5362986480000,
        bond_types::HC, 4.1025498300000,
        bond_types::CC, 4.7465390606595
    );

    // multiplicative constants by bond type in Angstroms
    constexpr double A = detail::get_bond_constant<bond_type>(
        bond_types::HH, 32.817355747000,
        bond_types::HC, 149.94098723000,
        bond_types::CC, 10953.544162170
    );

    // scaling constants for the inverse r potential term in eV
    constexpr double Q = detail::get_bond_constant<bond_type>(
        bond_types::HH, 0.3704714870450,
        bond_types::HC, 0.3407757280000,
        bond_types::CC, 0.3134602960833
    );

////////////////////////////////////////////////////////////////////////////////
// actual potential and potential derivative calculation                      //
////////////////////////////////////////////////////////////////////////////////

    // V^R(r) = (1 + Q / r) A exp(-alpha r)
    return std::make_tuple(
        // V^R(r)
        A * detail::lk_exp(-alpha * r) * (1 + Q / r),
        // d/dr V^R(r)
        A * detail::lk_exp(-alpha * r) * (-Q / (r * r) - alpha * (1 + Q / r))
    );
}

namespace detail {

/// Coefficients for the G_C(theta) and G_H(theta) angular functions
constexpr double coeff_gch[7][6] = {
    // G_C(theta) coefficients
    {0.375450000000003, 1.40677647515210,  2.25437749446242,  2.031282890266615, 1.429711740681553, 0.502401399437291}, // G_C(theta),   [1:cos(0.6082 * pi)]
    {0.271856000000000, 0.48891648922574,  -0.43308173807375, -0.55967720823125, 1.272040452823824, -0.04005399574456}, // gamma(theta), [1:cos(0.6082 * pi)]
    {0.707277245734054, 5.67743584848898,  24.09701877749944, 57.59183195998685, 71.88287000287255, 36.27886067346253}, // G_C(theta),   [cos(0.6082 * pi):-0.5]
    {0.002599999999998, -1.09800000000001, -4.34600000000003, -6.83000000000002, -4.92800000000002, -1.34240000000001}, // G_C(theta),   [-0.5:-1]
    // G_H(theta) coefficients
    {19.06502493209876, 2.01775628395062,  -2.56642191975309, 3.291332234567896, -2.65356150617283, 0.837669975308644}, // G_H(theta),   [-0.5:1]
    {16.95344062500006, -21.0823875000010, -102.468300000002, -210.643230000034, -229.847130000024, -94.9946400000077}, // G_H(theta),   [-5/6:-0.5]
    {270.4568000100967, 1549.63580005608,  3781.771900123691, 4582.154400136413, 2721.430800074573, 630.6336000163169}  // G_H(theta),   [-1:-5/6]
};

/// calculate the quintic polynomial and its derivative with coefficients
/// corresponding to the given id
constexpr auto calc_gch(const std::size_t coeff_id,
    const double x)
{
    decltype(auto) cg = coeff_gch[coeff_id];

    double power = 1,
        value = cg[0],
        derivative = 0;

    for (int i = 1; i < 6; ++i)
    {
        derivative += i * power * cg[i];

        power *= x;
        value += power * cg[i];
    }

    return std::make_tuple(value, derivative);
}

} // namespace detail

/// carbon G(theta)
constexpr auto gtheta_c(const double cos_theta, const double N_t)
{
    // default case, cos(theta) between [-1:-0.5]
    std::size_t coeff_id = 3;

    if (cos_theta > -0.333412 /* [cos(0.6082 * pi):1] */)
    {
        if (N_t < 3.2)
            coeff_id = 1;
        else if (N_t > 3.7)
            coeff_id = 0;
        else
        {
            // switches between G and gamma based on number of connected atoms
            const auto G = detail::calc_gch(0, cos_theta),
                gamma = detail::calc_gch(1, cos_theta),
                // Q switches from [3.2:3.7]
                Q = switching_fn_nocheck(2 * (N_t - 3.2));

            return std::make_tuple(
                // value, switches from G to gamma as N_t goes from 3.2 -> 3.7
                val(G) + val(Q) * (val(gamma) - val(G)),
                // d/d(cos_theta)
                deriv(G) + val(Q) * (deriv(gamma) + deriv(G)),
                // d/d(N_t)
                2 * deriv(Q) * (val(gamma) - val(G))
            );
        }
    }
    else if (cos_theta > -0.5 /* [-0.5:cos(0.6082 * pi)] */)
        coeff_id = 2;

    // need to pad output due to no second derivative term
    const auto result = detail::calc_gch(coeff_id, cos_theta);
    return std::make_tuple(val(result), deriv(result), 0.0);
}

/// hydrogen G(theta)
constexpr auto gtheta_h(const double cos_theta)
{
    if (cos_theta > -0.5 /* [-0.5:1] */)
        return detail::calc_gch(4, cos_theta);
    else if (cos_theta > -5.0 / 6.0 /* [-0.5:-5/6] */)
        return detail::calc_gch(5, cos_theta);

    // [-1:-5/6]
    return detail::calc_gch(6, cos_theta);
}

} // namespace test

#endif // SP2_TEST_HPP_HPP
