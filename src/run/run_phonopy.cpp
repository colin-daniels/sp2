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

#include <fstream>
#include <sstream>
#include <iomanip>
#include <phonopy/phonopy_io.hpp>

using namespace std;
using namespace sp2;

void relax_structure(structure_t &structure, run_settings_t rset);

int generate_displacements(phonopy::phonopy_settings_t pset);
int generate_force_sets(run_settings_t rset);

int write_irreps(phonopy::phonopy_settings_t pset);
int generate_bands(phonopy::phonopy_settings_t pset);
int generate_eigs(phonopy::phonopy_settings_t pset);
int generate_dos(phonopy::phonopy_settings_t pset);
/// write xyz animation file anim.xyz for the specified band_id (1 indexed)
int write_anim(phonopy::phonopy_settings_t pset, int mode_id);

vector<pair<double, vector<vec3_t>>> read_eigs();



void plot_modes(string filename, airebo::system_control_t &sys,
    vector<pair<double, vector<vec3_t>>> modes);

int write_spectra(sp2::phonopy::phonopy_settings_t pset,
    vector<pair<double, vector<vec3_t>>> modes, structure_t structure,
    const string &filename);

int write_raman_active_anim(sp2::phonopy::phonopy_settings_t pset,
    const std::vector<std::pair<double, double>> &spectra);

void write_log(string filename, string desc);

void output_xml(string filename, const vector<double> &gradient);


int sp2::run_phonopy(const run_settings_t &settings, MPI_Comm)
{
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
        if (generate_force_sets(settings) != 0)
        {
            std::cout << "Failed to generate force sets using phonopy."
                      << std::endl;
            return EXIT_FAILURE;
        }
    }

    write_log(settings.log_filename, "Computing Force Constants");
    if (settings.phonopy_settings.calc_raman)
    {
        if (generate_eigs(settings.phonopy_settings) != 0 ||
            write_irreps(settings.phonopy_settings) != 0)
            return EXIT_FAILURE;

        // read eigenmodes/etc
        auto modes = read_eigs();

        cout << "num modes: " << modes.size() << endl;
        if (!modes.empty())
            cout << "num eigs: " << modes[0].second.size() << endl;

        airebo::system_control_t sys;
        sys.init(structure);
        plot_modes("modes.dat", sys, modes);

        write_spectra(settings.phonopy_settings, modes, structure,
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


void output_xml(string filename, const vector<double> &gradient)
{
    ofstream outfile(filename);
    if (!outfile)
        return;

    outfile.precision(15);
    outfile << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"
               "<modeling>\n"
               " <generator>\n"
               "  <i name=\"program\" type=\"string\">vasp</i>\n"
               "  <i name=\"version\" type=\"string\">5.4.1</i>\n"
               " </generator>\n"
               " <calculation>\n"
               "  <varray name=\"forces\">\n";

    // we want forces, so -gradient
    for (auto v : sp2::dtov3(gradient))
        outfile << "<v>" << -v.x << ' ' << -v.y << ' ' << -v.z << "</v>\n";

    outfile << "  </varray>\n"
               " </calculation>\n"
               "</modeling>\n";
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

    minimize::acgsd_settings_t settings;
    settings.gradient_tolerance = 1e-5;

    auto minimize = [&](auto &&sys) {
        minimize::acgsd(sys.get_diff_fn(), sys.get_position(), settings);
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

    return system("phonopy -d disp.conf");
}

/// read all input poscar files, return vector of output filenames
vector<string> process_displacements(run_settings_t rset)
{
    auto get_ext = [](int i) {
        stringstream ss;
        ss << setw(3) << setfill('0') << i;
        return ss.str();
    };

    vector<string> output_filenames;

    std::unique_ptr<lammps::system_control_t> sys_lammps;
    std::unique_ptr<airebo::system_control_t> sys_rebo;

    std::function<std::vector<double>(void)> get_gradient;
    for (int i = 1; i < std::numeric_limits<int>::max(); ++i)
    {
        string input_file = "POSCAR-" + get_ext(i),
            output_file = "vasprun.xml-" + get_ext(i);

        structure_t structure;
        if (!io::file_exists(input_file) ||
            !io::read_structure(input_file, structure, file_type::POSCAR))
        {
            break;
        }

        if (!get_gradient)
        {
            auto grad_fn = [&](auto &&sys) -> std::vector<double> {
                sys.set_position(structure.positions);
                sys.update();
                return sys.get_gradient();
            };

            switch (rset.potential)
            {
            case potential_type::LAMMPS:
                sys_lammps = std::make_unique<lammps::system_control_t>(
                    structure, rset.lammps_settings);

                get_gradient = [&]{return grad_fn(*sys_lammps);};
                break;
            case potential_type::REBO:
                sys_rebo = std::make_unique<airebo::system_control_t>(
                    structure);

                get_gradient = [&]{return grad_fn(*sys_rebo);};
            default:
                return {};
            }
        }

        output_xml(output_file, get_gradient());
        output_filenames.push_back(output_file);
    }

    return output_filenames;
}

int generate_force_sets(run_settings_t rset)
{
    // process displacements
    auto force_filenames = process_displacements(rset);

    // generate force sets, "phonopy -f file1 file2 ..."
    string command = "phonopy -f";
    for (string file : force_filenames)
        command += " " + file;

    return system(command.c_str());
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
            << "FORCE_CONSTANTS = "
                << (io::file_exists("FORCE_CONSTANTS") ? "READ\n" : "WRITE\n");

    outfile.close();

    outfile.open("band_labels.txt");
    for (auto qp : pset.qpoints)
        outfile << '\t' << qp.label;

    outfile.close();

    return system("phonopy band.conf");
}

int write_irreps(phonopy::phonopy_settings_t pset)
{
    ofstream outfile("irreps.conf");
    outfile << "DIM = " << pset.supercell_dim[0]
            << ' ' << pset.supercell_dim[1]
            << ' ' << pset.supercell_dim[2] << '\n'
            << "IRREPS = 0 0 0 1e-3\n"
            << "FORCE_CONSTANTS = "
                << (io::file_exists("FORCE_CONSTANTS") ? "READ\n" : "WRITE\n");

    outfile.close();

    return system("phonopy irreps.conf");
}

int generate_eigs(phonopy::phonopy_settings_t pset)
{
    ofstream outfile("eigs.conf");

    outfile << "DIM = " << pset.supercell_dim[0]
                << ' ' << pset.supercell_dim[1]
                << ' ' << pset.supercell_dim[2] << '\n'
            << "BAND = 0 0 0   1/2 0 0\n"
            << "BAND_POINTS = 1\n"
            << "FORCE_CONSTANTS = "
                << (io::file_exists("FORCE_CONSTANTS") ? "READ\n" : "WRITE\n")
            << "EIGENVECTORS = .TRUE.\n";


    outfile.close();

    return system("phonopy eigs.conf");
}

int generate_dos(phonopy::phonopy_settings_t pset)
{
    ofstream outfile("dos.conf");

    outfile << "DIM = " << pset.supercell_dim[0]
                << ' ' << pset.supercell_dim[1]
                << ' ' << pset.supercell_dim[2] << '\n'
            << "MESH = 7 7 7\n"
            << "DOS = .TRUE.\n"
            << "DOS_RANGE = 0 90 0.1\n";

    outfile.close();

    return system("phonopy dos.conf");
}

int write_anim(phonopy::phonopy_settings_t pset, int mode_id)
{
    ofstream outfile("anim.conf");

    outfile << "DIM = " << pset.supercell_dim[0]
            << ' ' << pset.supercell_dim[1]
            << ' ' << pset.supercell_dim[2] << '\n'
            << "FORCE_CONSTANTS = "
            << (io::file_exists("FORCE_CONSTANTS") ? "READ\n" : "WRITE\n")
            << "ANIME_TYPE = XYZ\n"
            << "ANIME = " << mode_id << " 5 20\n";


    outfile.close();
    return system("phonopy anim.conf");
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
        for (int i = 0; i < mode.second.size(); ++i)
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


int write_spectra(sp2::phonopy::phonopy_settings_t pset,
    vector<pair<double, vector<vec3_t>>> modes, structure_t structure,
    const string &filename)
{
    constexpr double temperature = 293.15; // room temp
#warning temporary
    auto gaussian = [](double x, double peak)
    {
        // full width at half maximum for gaussian broadening
        constexpr double fwhm = 14;

        constexpr double sigma = fwhm / std::sqrt(std::log(256)),
            prefactor = 1 / (sigma * std::sqrt(2 * M_PI)),
            exp_denom = (2 * sigma * sigma);

        return prefactor * std::exp(-(x - peak) * (x - peak) / exp_denom);
    };

    vec3_t pol_vecs[3] = {
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1}
    };

    // phonopy outputs non-mass-normalized eigenvectors, so we
    // need to normalize them
    for (auto &mode : modes)
    {
        auto na = structure.types.size();
        for (std::size_t i = 0; i < na; ++i)
            if (structure.types[i] == atom_type::CARBON)
                mode.second[i] /= std::sqrt(12);
    }

    std::vector<std::pair<double, double>> spectra;
    if (pset.calc_raman_backscatter_avg)
    {
        spectra = sp2::phonopy::raman_spectra_avg(true, temperature,
            modes, structure);
    }
    else
    {
        spectra = sp2::phonopy::raman_spectra(
            pol_vecs[pset.polarization_axes[0]],
            pol_vecs[pset.polarization_axes[1]],
            temperature, modes, structure
        );
    }

    double maxi = 0,
        maxf = 0;
    for (auto mode : spectra)
    {
        maxi = std::max(mode.second, maxi);
        maxf = std::max(mode.first,  maxf);
    }

    constexpr int n_bins = 1200;
    const double bin_max = 1.1 * maxf,
        bin_step = bin_max / n_bins;
    std::vector<double> bins(n_bins, 0);

    auto irrep_labels = sp2::phonopy::read_irreps();

    ofstream outfile(filename);
    for (unsigned int i = 0; i < spectra.size(); ++i)
    {
        auto &shift = spectra[i];
        shift.second /= maxi;

        outfile << shift.first << ' '
                << shift.second << ' '
                << shift.second * maxi << ' '
                << irrep_labels[i] << endl;

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

    // write animation files
    if (pset.write_raman_active_anim)
        return write_raman_active_anim(pset, spectra);

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

int write_raman_active_anim(sp2::phonopy::phonopy_settings_t pset,
    const std::vector<std::pair<double, double>> &spectra)
{
    std::vector<unsigned int> raman_active_ids;
    for (auto i = 0u; i < spectra.size(); ++i)
    {
        if (spectra[i].second < pset.raman_active_anim_cutoff)
            continue;

        // phonopy starts counting bands at 1
        int mode_id = i + 1;

        // phonopy counts bands starting at 1
        if (write_anim(pset, mode_id) != 0 || !io::file_exists("anime.xyz"))
        {
            std::cerr << "Phonopy failed to write animation file for mode id "
                      << mode_id << ", quitting." << std::endl;

            return EXIT_FAILURE;
        }

        io::copy_file("anime.xyz", "anim_" + std::to_string(mode_id) + ".xyz");
    }

    return EXIT_SUCCESS;
}
