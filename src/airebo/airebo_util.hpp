#ifndef AIREBO_UTIL_HPP_INCLUDED
#define AIREBO_UTIL_HPP_INCLUDED

/// \file airebo_util.hpp
/// \brief Utility header with AIREBO-specific \cite stuart2000reactive potential and force functions

#include <cmath>

namespace sp2
{
    namespace airebo
    {
        /// constant for the LJ potential in lj_12_6()
        static const double epsilon[4] = {0.0015, sqrt(0.00284 * 0.0015), sqrt(0.00284 * 0.0015), 0.00284};

        /// constant for the LJ potential in lj_12_6()
        static const double sigma[4] = {2.65, 3.025, 3.025, 3.40};
        /// pre-exponentiated constants for the LJ potential in lj_12_6()
        static const double sig_6[4] =   {pow(sigma[0], 6),  pow(sigma[1], 6),  pow(sigma[2], 6),  pow(sigma[3], 6)};
        /// pre-exponentiated constants for the LJ potential in lj_12_6()
        static const double sig_12[4] = {pow(sigma[0], 12), pow(sigma[1], 12), pow(sigma[2], 12), pow(sigma[3], 12)};

        /// minimum switching distances for the LJ potential
        static const double lj_D_min[4] = {sigma[0], sigma[1], sigma[2], sigma[3]};
        /// maximum switching distances for the LJ potential
        static const double lj_D_max[4] = {pow(2, 1/6.0) * sigma[0], pow(2, 1/6.0) * sigma[1], pow(2, 1/6.0) * sigma[2], pow(2, 1/6.0) * sigma[3]};

        /// minimum switching bond order values for the LJ potential
        static const double lj_b_min[4] = {0.32, 0.90, 0.90, 0.81};
        /// maximum switching bond order values for the LJ potential
        static const double lj_b_max[4] = {0.42, 0.75, 0.75, 0.77};

        static inline void lj_t_b(const double b_ij, const int b_type, double &val, double &deriv)
        {
            if (b_ij < lj_b_min[b_type])
            {
                val = 1;
                deriv = 0;
            }
            else if (b_ij > lj_b_max[b_type])
            {
                val = 0;
                deriv = 0;
            }
            else
            {
                const double t = (b_ij - lj_b_min[b_type]) / (lj_b_max[b_type] - lj_b_min[b_type]);
                val = 1 - t * t * (3 - 2 * t);
                deriv = 6 * t * (t - 1) / (lj_b_max[b_type] - lj_b_min[b_type]);
            }
        }

        /// Lennard-Jones potential for van der Waals forces in AIREBO
        /// \param bond_length_sq const double bond length squared
        /// \param b_type const int bond type id
        /// \param lj double& output, lennard jones potential value
        /// \param d_lj double& output, derivative in the bond direction
        static inline void lj_12_6(const double bond_length_sq, const int b_type, double &lj, double &d_lj)
        {
            const double r6 = bond_length_sq * bond_length_sq * bond_length_sq;

            lj = 4 * epsilon[b_type] * (sig_12[b_type] / (r6 * r6) - sig_6[b_type] / r6);
            d_lj = -24 * epsilon[b_type] * (2 * sig_12[b_type] / (r6 * r6) - sig_6[b_type] / r6) / bond_length_sq;
        }

        /// switching function for LJ distances, assumes inputs are in range
        static inline void lj_S_r(const double bond_length, const int b_type, double &val, double &deriv)
        {
//            deriv = 0;
//            if (bond_length < lj_D_min[b_type])
//                val = 1;
//            else if (bond_length > lj_D_max[b_type])
//                val = 0;
//            else
//            {
                const double t = (bond_length - lj_D_min[b_type]) / (lj_D_max[b_type] - lj_D_min[b_type]);
                val = 1 - t * t * (3 - 2 * t);
                deriv = 6 * t * (t - 1) / (lj_D_max[b_type] - lj_D_min[b_type]);
//            }
        }

        /// switching function for LJ bond order
        static inline void lj_S_b(const double bond_order, const int b_type, double &val, double &deriv)
        {
            deriv = 0;
            if (bond_order < lj_b_min[b_type])
                val = 1;
            else if (bond_order > lj_b_max[b_type])
                val = 0;
            else
            {
                const double t = (bond_order - lj_b_min[b_type]) / (lj_b_max[b_type] - lj_b_min[b_type]);
                val = 1 - t * t * (3 - 2 * t);
                deriv = 6 * t * (t - 1) / (lj_b_max[b_type] - lj_b_min[b_type]);
            }
        }

        /// epsilon for single bond torsion
        static const double single_dh_epsilon[4] = {0.1250, 0.1787, 0.1787, 0.3079};

//        /// single bond torsion
//        inline void single_dh(const int type_a, const int type_d)
//        {
//            // V_tors = single_dh_epsilon[get_btype(type_a, type_d)] * (256/405 * (Theta_ijkl / 2)^10 - 0.1)
//        }
    }
}

#endif // AIREBO_UTIL_HPP_INCLUDED