#ifndef SP2_AIREBO_UTILITY_FUNCTIONS_HPP
#define SP2_AIREBO_UTILITY_FUNCTIONS_HPP

/// \file rebo_util.hpp
/// \brief Utility header with bond potential and force calculation functions for scpp::rebo::system_control_t

#include "common/util/blas.hpp"

namespace sp2 {
namespace airebo {

namespace {
/// calculate the exponential lookup table
std::vector<double> calc_exp_table(int n_points, double range)
{
    std::vector<double> result(n_points, 0);
    for (int i = 0; i < n_points; ++i)
        result[i] = std::exp(((i + 0.5) / n_points) * range);

    return result;
}

const int n_lk = 8192 * 32;
const double range_lk = -4.7204523127 * 2.1;
const std::vector<double> exp_lk = calc_exp_table(n_lk, range_lk);
} // anonymous namespace

/// exponential lookup table function, valid in the range [-9.9:0]
__attribute__((always_inline)) static inline double lk_exp(double x)
{
    const double dst = x - (std::trunc(n_lk * (x / range_lk)) + 0.5) * range_lk / n_lk;
    return ((1.0 + dst) + dst * dst * 0.5) * exp_lk[static_cast<int>(n_lk * (x / range_lk))];
}

///////////////
// constants for use in bond functions
///////////////

/// PI
static const double pi = std::atan2(0, -1);

/// minimum cutoff distance for bond types
static const double D_min[4] = {1.1, 1.3, 1.3, 1.7};
/// maximum cutoff distance for bond types
static const double D_max[4] = {1.7, 1.8, 1.8, 2.0};
/// (max - min) cutoff distance for bond types
static const double D_range[4] = {
    D_max[0] - D_min[0],
    D_max[1] - D_min[1],
    D_max[2] - D_min[2],
    D_max[3] - D_min[3]
};

/// scaling constants for the attractive potential, eV
static const double B[4][3] = {
    { 29.6325930000, 0, 0},
    { 32.3551866587, 0, 0},
    { 32.3551866587, 0, 0},
    {12388.79197798, 17.56740646509, 30.71493208065}
};

/// exponential constants for the attractive potential, Angstroms^-1
static const double beta[4][3] = {
    {1.71589217000, 0, 0},
    {1.43445805925, 0, 0},
    {1.43445805925, 0, 0},
    {4.72045231270, 1.4332132499, 1.3826912506}
};

/// term for the repulsive potential
static const double alpha[4] = {3.536298648000,   4.10254983,   4.10254983, 4.7465390606595};  // Angstroms^-1
/// term for the repulsive potential
static const double A[4]     = {32.81735574700, 149.94098723, 149.94098723, 10953.544162170};  // Angstroms
/// term for the repulsive potential
static const double Q[4]     = {0.370471487045,  0.340775728,  0.340775728, 0.3134602960833};  // eV

/// for the G(theta) sum
static const double lambda_ijk[4][4] = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 4}};

///////////////
// bond functions
///////////////

/// interaction cutoff function
static inline void cutoff_fn(const double &bond_length,
    const int b_type, double &cutoff, double &d_cutoff)
{
    if (bond_length < D_min[b_type])
        cutoff = 1;
    else if (bond_length < D_max[b_type])
    {
        double x = (bond_length - D_min[b_type]) / D_range[b_type];
        cutoff = (1 + 2 * x) * (1 - x) * (1 - x);
        d_cutoff = 6 * x * (x - 1) / D_range[b_type];
    }
}

/// attractive and repulsive potential terms
static inline void ar_ptnl(const double &bond_length, const int b_type,
    double &rep_ptnl, double &d_rep_ptnl, double &attr_ptnl, double &d_attr_ptnl)
{
    // repulsive potential
    rep_ptnl    = A[b_type] * std::exp(-alpha[b_type] * bond_length) * (1 + Q[b_type] / bond_length);
    d_rep_ptnl  = -A[b_type] * std::exp(-alpha[b_type] * bond_length) * (Q[b_type] / (bond_length * bond_length)
                                                                    + alpha[b_type] * (1 + Q[b_type] / bond_length));

    // attractive potential
    if (b_type == static_cast<int>(bond_type::CARBON_CARBON))
    {
        attr_ptnl   =  B[b_type][0] * exp(-beta[b_type][0] * bond_length) +
                       B[b_type][1] * exp(-beta[b_type][1] * bond_length) +
                       B[b_type][2] * exp(-beta[b_type][2] * bond_length);
        d_attr_ptnl = -beta[b_type][0] * B[b_type][0] * exp(-beta[b_type][0] * bond_length) +
                      -beta[b_type][1] * B[b_type][1] * exp(-beta[b_type][1] * bond_length) +
                      -beta[b_type][2] * B[b_type][2] * exp(-beta[b_type][2] * bond_length);
    }
    else
    {
        attr_ptnl   =  B[b_type][0] * exp(-beta[b_type][0] * bond_length);
        d_attr_ptnl = -beta[b_type][0] * B[b_type][0] * exp(-beta[b_type][0] * bond_length);
    }
}

/// attractive and repulsive potential terms lookup table ver
static inline void lk_ar_ptnl(const double &bond_length, const int b_type, double &rep_ptnl,
    double &d_rep_ptnl, double &attr_ptnl, double &d_attr_ptnl)
{
    // repulsive potential
    rep_ptnl    = A[b_type] * lk_exp(-alpha[b_type] * bond_length) * (1 + Q[b_type] / bond_length);
    d_rep_ptnl  = -A[b_type] * lk_exp(-alpha[b_type] * bond_length) * (Q[b_type] / (bond_length * bond_length)
                                                                       + alpha[b_type] * (1 + Q[b_type] / bond_length));

    // attractive potential
    if (b_type == static_cast<int>(bond_type::CARBON_CARBON))
    {
        attr_ptnl   =  B[b_type][0] * lk_exp(-beta[b_type][0] * bond_length) +
                       B[b_type][1] * lk_exp(-beta[b_type][1] * bond_length) +
                       B[b_type][2] * lk_exp(-beta[b_type][2] * bond_length);
        d_attr_ptnl = -beta[b_type][0] * B[b_type][0] * lk_exp(-beta[b_type][0] * bond_length) +
                      -beta[b_type][1] * B[b_type][1] * lk_exp(-beta[b_type][1] * bond_length) +
                      -beta[b_type][2] * B[b_type][2] * lk_exp(-beta[b_type][2] * bond_length);
    }
    else
    {
        attr_ptnl   =  B[b_type][0] * lk_exp(-beta[b_type][0] * bond_length);
        d_attr_ptnl = -beta[b_type][0] * B[b_type][0] * lk_exp(-beta[b_type][0] * bond_length);
    }
}

/// conjugation switching function
static inline double F_conj(const double x_ik)
{
    if (x_ik < 2)
        return 1;
    else if (x_ik < 3)
        return (x_ik - 3) * (x_ik - 3) * (2 * x_ik - 3);

    return 0;
}

/// conjugation switching function
static inline double dF_conj(const double x_ik)
{
    if (x_ik > 2 && x_ik < 3)
        return 6 * (x_ik - 2) * (x_ik - 3);

    return 0;
}

/// G_C(theta) coefficients
static double coeff_gc[4][6] = {
    {0.3754500000000032, 1.406776475152106, 2.254377494462427, 2.031282890266615, 1.429711740681553, 0.5024013994372916},           // G(theta),        [1:cos(0.6082 * pi)]
    {0.2718560000000009, 0.4889164892257488, -0.4330817380737564, -0.5596772082312531, 1.272040452823824, -0.04005399574456492},    // gamma(theta),    [1:cos(0.6082 * pi)]
    {0.7072772457340548, 5.677435848488983, 24.09701877749944, 57.59183195998685, 71.88287000287255, 36.27886067346253},            // G(theta),        [cos(0.6082 * pi):-0.5]
    {0.002599999999998603, -1.09800000000001, -4.34600000000003, -6.830000000000027, -4.928000000000027, -1.34240000000001}         // G(theta),        [-0.5:-1]
};

/// G_H(theta) coefficients
static double coeff_gh[3][6] = {
    {19.06502493209876, 2.01775628395062, -2.56642191975309, 3.291332234567896, -2.653561506172839, 0.8376699753086442},        // G(theta),    [-0.5:1]
    {16.95344062500006, -21.08238750000106, -102.4683000000023, -210.6432300000345, -229.8471300000242, -94.99464000000779},    // G(theta),    [-5/6:-0.5]
    {270.4568000100967, 1549.635800056089, 3781.771900123691, 4582.154400136413, 2721.430800074573, 630.6336000163169}          // G(theta),    [-1:-5/6]
};

/// G_C(theta) function
static inline void gtheta_c(const double x, const double N_t, double &gtheta, double &dgtheta, double &dgtheta_dq)
{
    const double inp[6] = {1, x, x*x, (x*x)*x, (x*x)*(x*x), (x*x)*(x*x*x)},
        d_inp[6] = {0, 1, 2*x, 3*(x*x), 4*(x*x*x), 5*(x*x)*(x*x)};

    if (x > -0.3334119772160521)
    {
        if (N_t > 3.7)
        {
            gtheta = loop_ddot<6>(inp, coeff_gc[0]);
            dgtheta = loop_ddot<6>(d_inp, coeff_gc[0]);
        }
        else if (N_t < 3.2)
        {
            gtheta = loop_ddot<6>(inp, coeff_gc[1]);
            dgtheta = loop_ddot<6>(d_inp, coeff_gc[1]);
        }
        else
        {
            double G = loop_ddot<6>(inp, coeff_gc[0]),
                dG = loop_ddot<6>(d_inp, coeff_gc[0]);

            double gamma = loop_ddot<6>(inp, coeff_gc[1]),
                dgamma = loop_ddot<6>(d_inp, coeff_gc[1]);

            double val_Q = 0.5 * (1 + cos(2 * pi * (N_t - 3.2))),
                val_dQ = - pi * sin(2 * pi * (N_t - 3.2));

            gtheta = G + val_Q * (gamma - G);
            dgtheta = dG + val_Q * (dgamma - dG);
            dgtheta_dq = val_dQ * (gamma - G);
        }
    }
    else if (x > -0.5)
    {
        gtheta = loop_ddot<6>(inp, coeff_gc[2]);
        dgtheta = loop_ddot<6>(d_inp, coeff_gc[2]);
    }
    else
    {
        gtheta = loop_ddot<6>(inp, coeff_gc[3]);
        dgtheta = loop_ddot<6>(d_inp, coeff_gc[3]);
    }
}

/// q(N_t) function for G_C(theta)
static inline void gtheta_q(const double N_t, double &val_Q, double &val_dQ)
{
    if (N_t < 3.2)
        val_Q = 1;
    else if (N_t > 3.7)
        val_Q = 0;
    else
    {
        const double x = (N_t - 3.2) / 0.5;
        val_Q = 2*x*x*x - 3*x*x + 1;
        val_dQ = (3*x*x - 6*x) * 2;
    }
}

/// hydrogen G(theta)
static inline void gtheta_H(const double x, double &gtheta, double &dgtheta)
{
    const double inp[6] = {1, x, x*x, x*x*x, x*x*x*x, x*x*x*x*x},
        d_inp[6] = {0, 1, 2*x, 3*x*x, 4*x*x*x, 5*x*x*x*x};

    const int id = x > -0.5 ? 0 : (x > -5.0/6.0 ? 1 : 2);
    gtheta = loop_ddot<6>(inp, coeff_gh[id]);
    dgtheta = loop_ddot<6>(d_inp, coeff_gh[id]);
}

/// carbon G(theta)
static inline void gtheta_C(const double x, const double N_t, double &gtheta, double &dgtheta, double &dgtheta_dq)
{
    // reset Q derivative
    dgtheta_dq = 0;

    if (x < -0.3333132475682373964842285) // theta > 0.6082 pi
    {
        // normal g(theta)
        if (x > -0.5)
        {
            gtheta = 0.693995817040389305703 + 5.51054053822947051689 * x + 23.26271593568917777866 * x*x + 55.5289880598589187228 * x*x*x
                     + 69.3649928712774536316 * x*x*x*x + 35.0653785222535841135 * x*x*x*x*x;
            dgtheta = 5.51054053822947051689 + 46.5254318713783555573 * x + 166.586964179576756169 * x*x
                      + 277.459971485109814527 * x*x*x + 175.326892611267920568 * x*x*x*x;
        }
        else
        {
            gtheta = 0.0026000000000000000000000 - 1.0980000000000000000000 * x - 4.346000000000000000000 * x*x - 6.830000000000000000000 * x*x*x
                     - 4.928000000000000000000 * x*x*x*x - 1.3424000000000000000000 * x*x*x*x*x;
            dgtheta = -1.098000000000000000000 - 8.69200000000000000000 * x - 20.4900000000000000000 * x*x
                      - 19.7120000000000000000 * x*x*x - 6.71200000000000000000 * x*x*x*x;
        }

    }
    else
    {
        if (N_t > 3.7)
        {
            // normal g(theta)
            gtheta = 0.375540000000000000000 + 1.40782654560568457227 * x + 2.25523991247825179535 * x*x + 2.02495229272588921688 * x*x*x
                     + 1.42412555038930119782 * x*x*x*x + 0.512315698800873217672 * x*x*x*x*x;
            dgtheta = 1.40782654560568457227 + 4.51047982495650359070 * x + 6.07485687817766765065 * x*x
                      + 5.69650220155720479129 * x*x*x + 2.56157849400436608836 * x*x*x*x;
        }
        else if (N_t < 3.2)
        {
            // gamma(theta)
            gtheta = 0.271856000000000000000 + 0.48929385553825719405 * x - 0.432797362614327722763 * x*x - 0.56183887680241383278 * x*x*x
                     + 1.27087433563367764688 * x*x*x*x - 0.0373879517551932853909 * x*x*x*x*x;
            dgtheta = 0.48929385553825719405 - 0.86559472522865544552 * x - 1.68551663040724149833 * x*x
                      + 5.0834973425347105875 * x*x*x - 0.186939758775966426955 * x*x*x*x;
        }
        else
        {
            double G = 0.375540000000000000000 + 1.40782654560568457227 * x + 2.25523991247825179535 * x*x + 2.02495229272588921688 * x*x*x
                       + 1.42412555038930119782 * x*x*x*x + 0.512315698800873217672 * x*x*x*x*x,
                dG = 1.40782654560568457227 + 4.51047982495650359070 * x + 6.07485687817766765065 * x*x
                     + 5.69650220155720479129 * x*x*x + 2.56157849400436608836 * x*x*x*x;

            double gamma = 0.271856000000000000000 + 0.48929385553825719405 * x - 0.432797362614327722763 * x*x - 0.56183887680241383278 * x*x*x
                           + 1.27087433563367764688 * x*x*x*x - 0.0373879517551932853909 * x*x*x*x*x,
                dgamma = 0.48929385553825719405 - 0.86559472522865544552 * x - 1.68551663040724149833 * x*x
                         + 5.0834973425347105875 * x*x*x - 0.186939758775966426955 * x*x*x*x;
            double val_Q = 0.5 * (1 + cos(2 * pi * (N_t - 3.2))),
                val_dQ = - pi * sin(2 * pi * (N_t - 3.2));

            gtheta = G + val_Q * (gamma - G);
            dgtheta = dG + val_Q * (dgamma - dG);
            dgtheta_dq = val_dQ * (gamma - G);
        }
    }
}

} // namespace rebo
} // namespace sp2

#endif // SP2_AIREBO_UTILITY_FUNCTIONS_HPP
