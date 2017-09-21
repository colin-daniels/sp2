#include <common/graph/ud_graph_t.hpp>
#include <common/math/vec3_t.hpp>
#include <airebo/system_control_t.hpp>
#include <common/util/random.hpp>
#include <common/io/structure.hpp>
#include "bond_polarization.hpp"

#include "common/math/mat3x3_t.hpp"
#include "common/math/vec3_util.hpp"

double raman_prefactor(
    double frequency,
    double temperature = 0,
    double incident_freq = 0
)
{
    // (hbar / k_b) in [K] per [cm-1]
    constexpr double hk = 0.22898852319;

    double bose_occupation = 1;
    if (temperature != 0)
        bose_occupation += 1 / (std::exp(hk * frequency / temperature) - 1);

    double multiplier = 1 / frequency;
    if (incident_freq != 0)
        multiplier *= std::pow(frequency - incident_freq, 4);

    return bose_occupation * multiplier;
}

double sp2::phonopy::raman_intensity(double frequency, double temperature,
    vec3_t incident, vec3_t scattered, const std::vector<vec3_t> &eigs,
    const std::vector<double> &masses, const graph::ud_graph_t &bond_graph,
    const std::vector<vec3_t> &bonds, const std::vector<atom_types> &types,
    const std::unordered_map<bond_types, pol_constant_t> &pol_constants)
{
    mat3x3_t tensor = raman_tensor(eigs, masses, bond_graph, bonds,
        types, pol_constants);


    const double prefactor = raman_prefactor(frequency, temperature),
        sum = dot(incident, tensor * scattered);

    return prefactor * (sum * sum);
}

double sp2::phonopy::raman_intensity_avg(bool backscatter, double frequency,
    double temperature, const std::vector<vec3_t> &eigs,
    const std::vector<double> &masses, const graph::ud_graph_t &bond_graph,
    const std::vector<vec3_t> &bonds, const std::vector<atom_types> &types,
    const std::unordered_map<bond_types, pol_constant_t> &pol_constants)
{
    const mat3x3_t tensor = raman_tensor(eigs, masses, bond_graph,
        bonds,
        types, pol_constants);

    // there was probably an easier way to do this, or a simple proof, given
    // the extremely simple answer
    //
    // random unit vectors in 3D can be generated by generating gaussian x, y, z
    //     v = (x, y, z) / sqrt(x^2 + y^2 + z^2)
    //       = (cos(phi) sin(theta), sin(phi) sin(theta), cos(theta))
    //
    // since we want the average of v^T (polarization tensor) v, we need to
    // find the expectation values for the matrix
    //
    //     1       (x1x2 a + x1y2 b + x1z2 c +
    // --------- *  y1x2 d + y1y2 e + y1z2 f +
    // (r1 r2)^2    z1x2 g + z1y2 h + z1z2 i)^2
    //
    // which ends up looking like the integral of
    //  = (elem)  e^(-(r1^2 + r2^2)/2) / (r1 r2)^2
    //  = (1/sqrt(2 pi))^6 * (elem / (r1 r2)^2) *  e^(-(r1^2 + r2^2)/2)
    //  = (1/sqrt(2 pi))^6 * g(theta1, phi1, theta2, phi2) * f(r1, r2)
    //
    // using [integral 0 to inf of (r^2 e^(-r^2/2) dr) = sqrt(pi / 2)] what
    // we have left is the integral of
    //  = 1 / (16 pi^2) (elem / (r1 r2)^2)
    //
    // so essentially all we need to integrate is
    //    [cos^2(phi) sin^2(theta), xx
    //     sin^2(phi) sin^2(theta), yy
    //     cos(phi) sin(phi) sin^2(theta), xy
    //     cos^2(theta), zz
    //     sin(phi) cos(theta) sin(theta), zy
    //     cos(phi) cos(theta) sin(theta)] zx
    //      * sin(theta) dtheta dphi
    //
    //     xx = 4 pi / 3
    //     yy = 4 pi / 3
    //     zz = 4 pi / 3
    //
    //     yx = 0
    //     yz = 0
    //     xz = 0
    //
    //     = [4/3 pi, 4/3 pi, 4/3 pi]
    //
    // back to the original equation we get
    //  = (4/3 pi)^2 / (16 pi^2)
    //  = 1 / 9
    //
    // for 2D (backscattering) its just:
    //  = 1 / (4 pi^2) (pi^2)
    //  = 1 / 4
    //
    const double prefactor = raman_prefactor(frequency, temperature);

    double sum = 0;

    const int lim = backscatter ? 2 : 3;
    for (int i = 0; i < lim; ++i)
        for (int j = 0; j < lim; ++j)
            sum += tensor[i][j] * tensor[i][j];

    if (backscatter)
        sum /= 4;
    else
        sum /= 9;

    return prefactor * sum;
}

sp2::mat3x3_t sp2::phonopy::raman_tensor(const std::vector<vec3_t> &eigs,
    const std::vector<double> &masses, const graph::ud_graph_t &bond_graph,
    const std::vector<vec3_t> &bonds, const std::vector<atom_types> &types,
    const std::unordered_map<bond_types, pol_constant_t> &pol_constants)
{
    // kronecker delta value
    constexpr double kdelta[3][3] = {
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1}
    };

    mat3x3_t tensor{};
    for (auto edge : bond_graph.edges())
    {
        const int bond_id = edge.id;
        const auto bond_type = btype(types[edge.a], types[edge.b]);

        // phonon eigenvector for this atom, need to mass normalize
        const vec3_t eig = eigs[edge.a] / std::sqrt(masses[edge.a]);

        // unit bond vector and length, used later
        const double len = bonds[bond_id].mag();
        const vec3_t rhat = bonds[bond_id] / len;

        // ignore bonds which have no corresponding polarization constants
        // specified in the input map
        if (!pol_constants.count(bond_type))
            throw std::runtime_error("No polarization constants specified "
                "for bond type " + enum_to_str(types[edge.a])
                                 + enum_to_str(types[edge.b]) + " "
                "between atoms " + std::to_string(edge.a) + " and "
                                 + std::to_string(edge.b));

        const auto& pc = pol_constants.at(bond_type);
        // check if bond is actually valid (via length)
        if (len > pc.max_len)
            continue;

        double const_one = pc.c1, // a || -   a |-
              dconst_one = pc.c2, // a'|| -   a'|-
              dconst_two = pc.c3; // a'|| + 2 a'|-

        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                tensor[i][j] -= dot(rhat, eig) * (
                      (dconst_two / 3) * kdelta[i][j]
                    + dconst_one * (rhat[i] * rhat[j] - kdelta[i][j] / 3)
                ) + (const_one / len) * (
                      (rhat[i] * eig[j] + rhat[j] * eig[i])
                    - 2 * rhat[i] * rhat[j] * dot(rhat, eig)
                );
            }
        }
    }

    return tensor;
}

std::vector<std::pair<double, double>> raman_spectra_impl(
    bool avg, bool backscatter,
    sp2::vec3_t incident, sp2::vec3_t scattered, double temperature,
    const std::vector<std::pair<double, std::vector<sp2::vec3_t>>> &modes,
    const std::vector<double> &masses,
    const sp2::structure_t &structure
)
{
    sp2::airebo::system_control_t sys;
    sys.init(structure);
    sys.update();

    auto bond_graph = sys.get_bond_control().get_graph();
    auto bond_deltas = sp2::dtov3(sys.get_bond_control().get_bond_deltas());

    // convert from old to new types for now
    const auto &types = structure.types;

    std::vector<std::pair<double, double>> result;
    for (auto mode : modes)
    {
        double frequency = mode.first;
        const auto &eigs = mode.second;

        // calculate raman intensity for given the incident/scattered light
        // polarization directions and the current eigenmode
        double intensity;
        if (!avg)
        {
            intensity = sp2::phonopy::raman_intensity(frequency, temperature,
                incident, scattered, eigs,
                masses, bond_graph, bond_deltas, types);
        }
        else
        {
            intensity = sp2::phonopy::raman_intensity_avg(backscatter,
                frequency, temperature, eigs, masses,
                bond_graph, bond_deltas, types);
        }

        result.emplace_back(frequency, intensity);
    }

    // note that we don't normalize to unity, the user must
    return result;
}

std::vector<std::pair<double, double>> sp2::phonopy::raman_spectra(
    vec3_t incident, vec3_t scattered, double temperature,
    const std::vector<std::pair<double, std::vector<vec3_t>>> &modes,
    const std::vector<double> &masses,
    const structure_t &structure
)
{
    return raman_spectra_impl(false, false, incident, scattered,
        temperature, modes, masses, structure);
}


std::vector<std::pair<double, double>> sp2::phonopy::raman_spectra_avg(
    bool backscatter, double temperature,
    const std::vector<std::pair<double, std::vector<vec3_t>>> &modes,
    const std::vector<double> &masses,
    const structure_t &structure
)
{
    return raman_spectra_impl(true, backscatter, vec3_t{}, vec3_t{},
        temperature, modes, masses, structure);
}
