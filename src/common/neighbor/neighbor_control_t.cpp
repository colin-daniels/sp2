#include "neighbor_control_t.hpp"
#include <gtest/gtest.h>
#include <common/math/mat3x3_t.hpp>
#include <common/util/random.hpp>

sp2::mat3x3_t gen_lattice(sp2::util::rng_t &rng, double r_min = 1, double r_max = 100)
{
    auto gen_vec = [&]{
        return sp2::random_vec3(rng.get_gen()) * rng.rand(r_min, r_max);
    };

    sp2::vec3_t lattice_vecs[3] = {
        gen_vec(), gen_vec(), gen_vec()
    };

    sp2::mat3x3_t lattice;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            lattice[j][i] = lattice_vecs[i][j];

    return lattice;
}

TEST(lattice, all)
{
    auto rng = []{
        sp2::util::rng_t rng;
        rng.seed_random();
        return rng;
    }();

    constexpr int nt = 0,
        n_points = 1e2;

    std::vector<sp2::vec3_t> points(n_points);
    for (int t = 0; t < nt; ++t)
    {
        auto mat = gen_lattice(rng);
        for (auto &p : points)
            p = sp2::vec3_t(rng.rand(0.0, 1.0), rng.rand(0.0, 1.0), rng.rand(0.0, 1.0));

        for (auto &a : points)
        {
            for (auto &b : points)
            {
                auto delta = b - a;
                for (auto &coord : delta)
                    coord -= std::round(coord);

                double round_dist = (mat * delta).mag_sq();

                double brute_force_dist = std::numeric_limits<double>::max();
                for (int i = -1; i <= 1; ++i)
                    for (int j = -1; j <= 1; ++j)
                        for (int k = -1; k <= 1; ++k)
                            brute_force_dist = std::min(brute_force_dist,
                                (mat * (delta + sp2::vec3_t(i, j, k))).mag_sq());

                ASSERT_DOUBLE_EQ(round_dist, brute_force_dist);
            }
        }

        std::cout << "t: " << t << std::endl;
    }
}
