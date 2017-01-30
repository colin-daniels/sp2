#include "symm/util.hpp"
#include "common/util/random.hpp"

#include <vector>
#include <cmath>
#include <fstream>

using namespace std;
using namespace sp2;

void gen_sol(double r1, double r2, double &vn, double &vp, const double alpha = 0)
{
    const double pi = atan2(0, -1);

    double cr1, cr2, sr1, sr2;
    sincos(2 * pi * r1, &sr1, &cr1);
    sincos(2 * pi * r2, &sr2, &cr2);

    vn = atan((cr2 - sqrt(-alpha*alpha + 2*alpha*cr1*sr2 - cr1*cr1*sr2*sr2 + sr1*sr1 + cr2*cr2)) / (alpha - cr1*sr2 + sr1)) / pi;
    vp = atan((cr2 + sqrt(-alpha*alpha + 2*alpha*cr1*sr2 - cr1*cr1*sr2*sr2 + sr1*sr1 + cr2*cr2)) / (alpha - cr1*sr2 + sr1)) / pi;

    vp = vp < 0 ? vp + 1 : vp;
    vn = vn < 0 ? vn + 1 : vn;
}

void symm::test()
{
    // const int n = 6;

    vector<int> delta_f(13);
    delta_f[6] = 0;

    delta_f[4] = 2310;
    delta_f[5] = 924;
    delta_f[7] = -660;
    delta_f[8] = -1155;
    delta_f[9] = -1540;
    delta_f[10] = -1848;
    delta_f[11] = -2100;
    delta_f[12] = -2310;

    // // octagons
    // for (int n = 2; n < 10; ++n)
    // {
    //     vector<int> path(n, 6);
    // }

    util::rng_t rng;
    ofstream outfile("surf.dat");
    for (int i = 0; i < 1000; ++i)
    {
        double x = rng.rand(),
            y = rng.rand(),
            z1, z2;

        gen_sol(y, x, z1, z2);
        outfile << x << '\t' << y << '\t' << z1 << '\n'
                << x << '\t' << y << '\t' << z2 << '\n';
    }
}


double symm::gyroid_ptnl(const double *pos, double *deriv, const double c)
{
    static const double pi = atan2(0, -1);
    double cx, cy, cz,
        sx, sy, sz;
    sincos(2 * pi * pos[0], &sx, &cx);
    sincos(2 * pi * pos[1], &sy, &cy);
    sincos(2 * pi * pos[2], &sz, &cz);

    double val = cx*sz + cz*sy + cy*sx - c;

    deriv[0] = 2 * val * (-sx * sz + cy * cx) * 2 * pi;
    deriv[1] = 2 * val * (-sy * sx + cz * cy) * 2 * pi;
    deriv[2] = 2 * val * (-sz * sy + cz * cx) * 2 * pi;

    return val * val;
}

// structure_t symm::get_irreducible(const structure_t &input,
//     const std::map<std::string, space_group_t> &groups)
// {
//     structure_t fundamental;

//     using ref_type = std::reference_wrapper<const space_group_t>;

//     vector<ref_type> possible;

//     for (auto &elem : groups)
//         possible.push_back(std::cref(elem.second));

//     return fundamental;
// }

// struct plane3_t
// {
//     vec3_t position,
//         normal;
// };

// vector<plane3_t> make_voronoi(const vector<vec3_t> &points, size_t point_idx)
// {

// }

// void test_irreducible()
// {
//     auto group = symm::read_groups("space_groups.txt")["I a -3 d"];

//     // starting lattice point
//     auto point = vec3_t(0, 0, 0);

//     // apply all symmetries
//     auto points = group.apply_symm(vector<vec3_t>{point});

//     // wrap points into the unit cell (components between 0 and 1)
//     for (auto &p : points)
//         for (auto &component : p)
//             component -= floor(component);

//     // translate entire set of points to get neighboring cells as well
//     // (only those with negative translations though)
//     vector<vec3_t> all_points;
//     auto add_offset = [&](int dx, int dy, int dz)
//     {
//         auto temp = points;
//         for (auto &p : temp)
//             all_points.push_back(p + vec3_t(dx, dy, dz));
//     };

//     for (int x = -1; x <= 0; ++x)
//         for (int y = -1; y <= 0; ++y)
//             for (int z = -1; z <= 0; ++z)
//                 add_offset(x, y, z);

// }
