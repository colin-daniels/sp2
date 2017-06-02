#include "run/run_types.hpp"
#include "common/io/util.hpp"
#include "common/io/structure.hpp"
#include "airebo/system_control_t.hpp"
#include "common/math/blas.hpp"
#include "common/minimize/minimize.hpp"
#include "common/math/vec3_t.hpp"
#include "common/math/vec3_util.hpp"
#include "common/json/json.hpp"
#include "phonopy/bond_polarization.hpp"
#include "common/util/modeling.hpp"

#ifdef SP2_ENABLE_LAMMPS
#include "lammps/lammps_interface.hpp"
#endif

#include <thread>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <phonopy/phonopy_io.hpp>
#include <common/math/rotations.hpp>

using namespace std;
using namespace sp2;

inline std::string get_force_constant_rw()
{
    return std::string("FORCE_CONSTANTS = ") + (
        io::file_exists("FORCE_CONSTANTS")      ? "READ\n" :
        io::file_exists("force_constants.hdf5") ? "READ\n" :
             "WRITE\n"
    );
}

inline int run_phonopy_with_args(const phonopy::phonopy_settings_t pset,
    std::string args)
{
    std::string command = "phonopy "
          "--tolerance=" + std::to_string(pset.symmetry_tol) + " "
        + "--hdf5 "
        + args;

    return system(command.c_str());
}

void relax_structure(structure_t &structure, run_settings_t rset);

int generate_displacements(phonopy::phonopy_settings_t pset);
void generate_force_sets(run_settings_t rset);

int write_irreps(phonopy::phonopy_settings_t pset);
int generate_bands(phonopy::phonopy_settings_t pset);
int generate_eigs(phonopy::phonopy_settings_t pset);
int generate_dos(phonopy::phonopy_settings_t pset);
/// write xyz animation file anim.xyz for the specified band_id (1 indexed)
int write_anim(phonopy::phonopy_settings_t pset, int mode_id);

vector<pair<double, vector<vec3_t>>> read_eigs();


void plot_modes(string filename, airebo::system_control_t &sys,
    vector<pair<double, vector<vec3_t>>> modes);

int write_spectra(run_settings_t rset,
    vector<pair<double, vector<vec3_t>>> modes, structure_t structure,
    const string &filename);

void write_log(string filename, string desc);

void remove_hydrogen(structure_t &structure,
    vector<pair<double, vector<vec3_t>>> *modes = nullptr)
{
    auto &pos = structure.positions;
    for (std::size_t i = 0; i < structure.types.size(); ++i)
    {
        if (structure.types[i] != atom_type::HYDROGEN)
            continue;

        pos.erase(pos.begin() + i);
        structure.types.erase(structure.types.begin() + i);
        if (modes)
        {
            for (auto &mode : *modes)
                mode.second.erase(mode.second.begin() + i);
        }

        i -= 1;
    }
}

std::vector<double> get_phonopy_masses(const sp2::structure_t &structure)
{
    std::vector<double> masses;
    masses.reserve(structure.types.size());

    for (const auto &type : structure.types)
    {
        switch (type)
        {
        case atom_type::CARBON:
            masses.push_back(12.0107);
            break;
        case atom_type::HYDROGEN:
            masses.push_back(1.00794);
            break;
        }
    }

    return masses;
}

std::string get_phonopy_mass_command(const phonopy::phonopy_settings_t &pset)
{
    if (pset.masses.empty())
        return "";

    std::string mass_setting = "MASS =";
    for (auto mass : pset.masses)
        mass_setting += " " + std::to_string(mass);

    return mass_setting + "\n";
}

int sp2::run_phonopy(const run_settings_t &settings_in, MPI_Comm)
{
    // copy settings so we can modify it
    auto settings = settings_in;
    // copy structure for convenience
    auto structure = settings.structure;

#ifdef SP2_DEBUG
    io::write_structure("input.xyz", structure);
#endif // SP2_DEBUG

    // relax
    write_log(settings.log_filename, "Relaxing Structure");
    relax_structure(structure, settings);

    sort_structure_types(structure);
    io::write_structure("POSCAR", structure, false, file_type::POSCAR);

#ifdef SP2_DEBUG
    io::write_structure("relaxed.xyz", structure);
#endif // SP2_DEBUG

    // set masses if none are present and raman is going to be calculated
    if (settings.phonopy_settings.masses.empty() &&
        settings.phonopy_settings.calc_raman)
    {
        settings.phonopy_settings.masses = get_phonopy_masses(structure);
    }

    if (settings.phonopy_settings.calc_displacements)
    {
        write_log(settings.log_filename, "Generating Displacements");
        if (generate_displacements(settings.phonopy_settings) != 0)
        {
            std::cout << "Failed to generate displacements using phonopy."
                      << std::endl;
            return EXIT_FAILURE;
        }
    }

    if (settings.phonopy_settings.calc_force_sets)
    {
        write_log(settings.log_filename, "Generating Force Sets");
        generate_force_sets(settings);
    }

    write_log(settings.log_filename, "Computing Force Constants");
    if (settings.phonopy_settings.calc_raman)
    {
        if (generate_eigs(settings.phonopy_settings) != 0)
        {
            std::cerr << "Failed to calculate eigenvectors/eigenvalues.\n";
            return EXIT_FAILURE;
        }

        if (settings.phonopy_settings.calc_irreps &&
            write_irreps(settings.phonopy_settings) != 0)
        {
            std::cerr << "Failed to calculate irreps.\n";
            return EXIT_FAILURE;
        }

        // read eigenmodes/etc
        auto modes = read_eigs();

        cout << "num modes: " << modes.size() << endl;
        if (!modes.empty())
            cout << "num eigs: " << modes[0].second.size() << endl;

        airebo::system_control_t sys;
        sys.init(structure);
        plot_modes("modes.dat", sys, modes);

        write_spectra(settings, modes, structure,
            "spectra.dat");
    }

    if (settings.phonopy_settings.calc_bands &&
        generate_bands(settings.phonopy_settings) != 0)
    {
        std::cout << "Failed to generate phonon band structure using phonopy."
                  << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

void relax_structure(structure_t &structure, run_settings_t rset)
{
    // construct system + minimize with default parameters
    if (rset.add_hydrogen)
    {
        airebo::system_control_t sys;
        sys.init(structure);
        sys.add_hydrogen();
        structure = sys.get_structure();
    }

    auto supercell = util::construct_supercell(structure,
        rset.phonopy_settings.supercell_dim[0],
        rset.phonopy_settings.supercell_dim[1],
        rset.phonopy_settings.supercell_dim[2]);

    auto minimize = [&](auto &&sys) {
        minimize::acgsd(sys.get_diff_fn(), sys.get_position(),
            rset.phonopy_settings.min_set);
        supercell = sys.get_structure();
    };

    switch (rset.potential)
    {
    case potential_type::LAMMPS:
        minimize(lammps::system_control_t(supercell, rset.lammps_settings));
        break;
    case potential_type::REBO:
        minimize(airebo::system_control_t(supercell));
        break;
    default:
        return;
    }

    // just get the positions from the first primitive cell of the
    // supercell and pass it back out
    supercell.positions.resize(structure.positions.size());

    structure.positions = supercell.positions;
}

int generate_displacements(phonopy::phonopy_settings_t pset)
{
    ofstream outfile("disp.conf");

    outfile << "DIM = " << pset.supercell_dim[0]
                << ' ' << pset.supercell_dim[1]
                << ' ' << pset.supercell_dim[2] << '\n'
            << "DISPLACEMENT_DISTANCE = " << pset.displacement_distance << "\n"
            << "CREATE_DISPLACEMENTS = .TRUE.\n";

    outfile.close();

    return run_phonopy_with_args(pset, "disp.conf");
}

/// read disp.yaml and calculate forces to output to FORCE_SETS
void generate_force_sets(run_settings_t rset)
{
    // read in structure and displacements from phonopy's disp.yaml
    structure_t structure;
    std::vector<std::pair<int, sp2::vec3_t>> displacements;

    auto disp_data = phonopy::read_displacements("disp.yaml");
    std::tie(structure, displacements) = disp_data;

    // store original position
    auto orig_pos = structure.positions;


    // pointers (and function) for switching potentials for gradient calc
    std::unique_ptr<lammps::system_control_t> sys_lammps;
    std::unique_ptr<airebo::system_control_t> sys_rebo;
    std::function<std::vector<vec3_t>(const std::vector<double>&)> get_forces;


    // go through each displacement and calculate forces
    std::vector<std::vector<vec3_t>> forces;
    for (auto disp : displacements)
    {
        // copy original positions, and displace the specified atom
        auto temp_pos = orig_pos;
        temp_pos[disp.first] += disp.second;

        if (!get_forces)
        {
            auto force_fn = [&](auto &pos, auto &&sys) -> std::vector<vec3_t> {
                sys.set_position(pos);
                sys.update();

                // get gradient, negate it for forces
                auto temp = sys.get_gradient();
                vscal(-1.0, temp);

                return dtov3(temp);
            };

            switch (rset.potential)
            {
            case potential_type::LAMMPS:
                sys_lammps = std::make_unique<lammps::system_control_t>(
                    structure, rset.lammps_settings);

                get_forces = [&](auto& pos){return force_fn(pos, *sys_lammps);};
                break;
            case potential_type::REBO:
                sys_rebo = std::make_unique<airebo::system_control_t>(
                    structure);

                get_forces = [&](auto& pos){return force_fn(pos, *sys_rebo);};
                break;
            default:
                return;
            }
        }

        forces.emplace_back(
            get_forces(v3tod(temp_pos))
        );
    }

    // output FORCE_SETS for phonopy
    phonopy::write_force_sets("FORCE_SETS", disp_data, forces);
}

int generate_bands(phonopy::phonopy_settings_t pset)
{
    ofstream outfile("band.conf");

    outfile << "DIM = " << pset.supercell_dim[0]
                << ' ' << pset.supercell_dim[1]
                << ' ' << pset.supercell_dim[2] << '\n'
            << "BAND =";

    for (auto qp : pset.qpoints)
        outfile << "  " << qp.x << ' ' << qp.y << ' ' << qp.z;

    outfile << "\n"
            << "BAND_POINTS = " << pset.n_samples << "\n"
            << get_phonopy_mass_command(pset)
            << get_force_constant_rw();

    outfile.close();

    outfile.open("band_labels.txt");
    for (auto qp : pset.qpoints)
        outfile << '\t' << qp.label;

    outfile.close();

    return run_phonopy_with_args(pset, "band.conf");
}

int write_irreps(phonopy::phonopy_settings_t pset)
{
    ofstream outfile("irreps.conf");
    outfile << "DIM = " << pset.supercell_dim[0]
            << ' ' << pset.supercell_dim[1]
            << ' ' << pset.supercell_dim[2] << '\n'
            << "IRREPS = 0 0 0 1e-2\n"
            << get_phonopy_mass_command(pset)
            << get_force_constant_rw();

    outfile.close();

    return run_phonopy_with_args(pset, "irreps.conf");
}

int generate_eigs(phonopy::phonopy_settings_t pset)
{
    ofstream outfile("eigs.conf");

    outfile << "DIM = " << pset.supercell_dim[0]
                << ' ' << pset.supercell_dim[1]
                << ' ' << pset.supercell_dim[2] << '\n'
            << "BAND = 0 0 0   1/2 0 0\n"
            << "BAND_POINTS = 1\n"
            << get_phonopy_mass_command(pset)
            << get_force_constant_rw()
            << "EIGENVECTORS = .TRUE.\n";


    outfile.close();

    return run_phonopy_with_args(pset, "eigs.conf");
}

int generate_dos(phonopy::phonopy_settings_t pset)
{
    ofstream outfile("dos.conf");

    outfile << "DIM = " << pset.supercell_dim[0]
                << ' ' << pset.supercell_dim[1]
                << ' ' << pset.supercell_dim[2] << '\n'
            << get_force_constant_rw()
            << get_phonopy_mass_command(pset)
            << "MESH = 7 7 7\n"
            << "DOS = .TRUE.\n"
            << "DOS_RANGE = 0 90 0.1\n";

    outfile.close();

    return run_phonopy_with_args(pset, "dos.conf");
}

int write_anim(phonopy::phonopy_settings_t pset, int mode_id)
{
    ofstream outfile("anim.conf");

    outfile << "DIM = " << pset.supercell_dim[0]
            << ' ' << pset.supercell_dim[1]
            << ' ' << pset.supercell_dim[2] << '\n'
            << get_force_constant_rw()
            << get_phonopy_mass_command(pset)
            << "ANIME_TYPE = XYZ\n"
            << "ANIME = " << mode_id << " 5 20\n";


    outfile.close();
    return run_phonopy_with_args(pset, "anim.conf");
}

vector<pair<double, vector<vec3_t>>> read_eigs()
{
    ifstream infile("band.yaml");
    if (!infile)
    {
        cout << "Error opening band.yaml" << endl;
        return {};
    }

    vector<pair<double, vector<vec3_t>>> modes;
    // are we currently reading a eigenmode
    bool reading_mode = false;

    int idx = 0;
    for (string line; getline(infile, line); line.clear())
    {
        if (line.find("frequency") == string::npos)
        {
            if (!reading_mode)
                continue;

            auto pos = line.find("[");
            if (pos == string::npos)
                continue;

            line = line.substr(pos + 1, line.size() - pos);
            line = line.substr(0, line.find("]") - 1);

            double comp_real = 0,
                comp_img = 0;
            stringstream ss(line);
            ss >> comp_real;
            ss >> comp_img;

            auto &vec = modes.back().second;
            if (idx == 0)
                vec.emplace_back(comp_real, 0, 0);
            else
                vec.back()[idx] = comp_real;

            idx = (idx + 1) % 3;
        }
        else
        {
            double frequency = 0;
            auto pos = line.find(":");
            line = line.substr(pos + 1, line.size() - pos - 1);

            stringstream ss(line);
            ss >> frequency;

            // phonopy outputs frequency in terahertz, convert to [cm^-1]
            frequency *= 33.35641;

            modes.emplace_back(frequency, vector<vec3_t>{});
            reading_mode = true;
        }
    }

    return modes;
}

void plot_modes(string filename, airebo::system_control_t &sys,
    vector<pair<double, vector<vec3_t>>> modes)
{
    ofstream outfile(filename);

    auto positions = dtov3(sys.get_position());
    auto bonds = dtov3(sys.get_bond_control().get_bond_deltas());

    for (auto &mode : modes)
    {
        for (std::size_t i = 0; i < mode.second.size(); ++i)
        {
            auto pos1 = positions[i],
                pos2 = mode.second[i];


            for (int k = 0; k < 3; ++k)
                outfile << pos1[k] << ' ';

            for (int k = 0; k < 3; ++k)
                outfile << pos2[k] << ' ';

            outfile << endl;
        }

        outfile << "\n\n";
    }
}

void renomalize(pair<double, vector<vec3_t>> &mode)
{
    double sum = 0;
    for (auto e : mode.second)
        sum += dot(e, e);

    for (auto &e : mode.second)
        e /= sum;
}

int write_spectra(run_settings_t rset,
    vector<pair<double, vector<vec3_t>>> modes, structure_t structure,
    const string &filename)
{
    constexpr double temperature = 293.15; // room temp
    auto gaussian = [](double x, double peak)
    {
        // full width at half maximum for gaussian broadening
        constexpr double fwhm = 15;

        const double sigma = fwhm / std::sqrt(std::log(256)),
            prefactor = 1 / (sigma * std::sqrt(2 * M_PI)),
            exp_denom = (2 * sigma * sigma);

        return prefactor * std::exp(-(x - peak) * (x - peak) / exp_denom);
    };

    vec3_t pol_vecs[3] = {
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1}
    };

    std::vector<std::pair<double, double>> spectra;
    if (rset.phonopy_settings.calc_raman_backscatter_avg)
    {
        spectra = phonopy::raman_spectra_avg(true,
            temperature, modes, rset.phonopy_settings.masses, structure);
    }
    else
    {
        spectra = sp2::phonopy::raman_spectra(
            pol_vecs[rset.phonopy_settings.polarization_axes[0]],
            pol_vecs[rset.phonopy_settings.polarization_axes[1]],
            temperature, modes, rset.phonopy_settings.masses, structure
        );
    }

    // ignore acoustic modes
    for (std::size_t i = 0; i < spectra.size() && spectra[i].first < 1; ++i)
        spectra[i].second = 0;

    double maxi = 0,
        maxf = 0;
    for (auto mode : spectra)
    {
        maxi = std::max(mode.second, maxi);
        maxf = std::max(mode.first,  maxf);
    }

    const std::size_t n_bins = std::max<std::size_t>(maxf, 1000);
    const double bin_max = 1.1 * maxf,
        bin_step = bin_max / n_bins;
    std::vector<double> bins(n_bins, 0);

    std::vector<std::string> irrep_labels;
    if (!rset.phonopy_settings.calc_irreps)
        irrep_labels.resize(spectra.size(), "None");
    else
        irrep_labels = sp2::phonopy::read_irreps();

    if (spectra.size() != irrep_labels.size())
        throw std::runtime_error("Number of characters read from irreps.yaml ("
            + std::to_string(irrep_labels.size()) +
            ") does not match number of bands ("
            + std::to_string(spectra.size()) + ")");

    ofstream outfile(filename);
    for (unsigned int i = 0; i < spectra.size(); ++i)
    {
        auto &shift = spectra[i];
        shift.second /= maxi;

        outfile << shift.first << ' '
                << shift.second << ' '
                << shift.second * maxi << ' '
                << irrep_labels[i] << ' '
                << (i + 1) << '\n';

        for (std::size_t j = 0; j < bins.size(); ++j)
            bins[j] += gaussian(j * bin_step, shift.first) * shift.second;
    }
    outfile.close();

    outfile.open("gauss_" + filename);
    double max_bin = *std::max_element(bins.begin(), bins.end());

    for (std::size_t i = 0; i < bins.size(); ++i)
        outfile << (i * bin_step) << ' '
                << bins[i] / max_bin << std::endl;

    outfile.close();

    if (!(rset.phonopy_settings.write_all_mode_anim ||
          rset.phonopy_settings.write_all_mode_gplot))
        return EXIT_SUCCESS;

    // write animation files
    std::vector<unsigned int> raman_active_ids;
    auto get_id_suffix = [id_max = spectra.size()](int id) {
        const int n_digits = static_cast<int>(std::log10(id_max)) + 1;

        std::ostringstream oss;
        oss << std::setw(n_digits) << std::setfill('0') << id;
        return oss.str();
    };

    for (auto i = 0u; i < spectra.size(); ++i)
    {
        // phonopy starts counting bands at 1
        int mode_id = i + 1;

        // gnuplot file
        if (rset.phonopy_settings.write_all_mode_gplot)
        {
            phonopy::draw_normal_mode(
                "mode_" + get_id_suffix(mode_id) + ".gplot",
                structure, modes[i]);
        }

        // animation file
        if (rset.phonopy_settings.write_all_mode_anim)
        {
            if (write_anim(rset.phonopy_settings, mode_id) != 0 ||
                !io::file_exists("anime.xyz"))
            {
                std::cerr
                    << "Phonopy failed to write animation file for mode id "
                    << mode_id << ", quitting." << std::endl;

                return EXIT_FAILURE;
            }

            io::copy_file("anime.xyz",
                "anim_" + get_id_suffix(mode_id) + ".xyz");
        }
    }

    return EXIT_SUCCESS;
}

void write_log(string filename, string desc)
{
    if (filename.empty())
        return;

    stringstream ss;
    ss << "{\"status\":\"" << desc << "\"}";

    string string = ss.str();
    ofstream outfile(filename);
    outfile.write(string.c_str(), string.size());
    outfile.close();
}
