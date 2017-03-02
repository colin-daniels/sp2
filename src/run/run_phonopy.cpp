#include "run/run_types.hpp"
#include "common/io/util.hpp"
#include "common/io/structure.hpp"
#include "airebo/system_control_t.hpp"
#include "common/util/blas.hpp"
#include "common/minimize/minimize.hpp"
#include "common/vec3_t.hpp"
#include "common/json/json.hpp"
#include "phonopy/bond_polarization.hpp"
#include "common/util/modeling.hpp"

#ifdef SP2_ENABLE_LAMMPS
#include "lammps/lammps_interface.hpp"
#endif

#include <fstream>

using namespace std;
using namespace sp2;

void relax_structure(structure_t &structure, run_settings_t rset);

int generate_displacements(phonopy::phonopy_settings_t pset);
int generate_force_sets(run_settings_t rset);

int generate_bands(phonopy::phonopy_settings_t pset);
int generate_eigs(phonopy::phonopy_settings_t pset);
int generate_dos(phonopy::phonopy_settings_t pset);

vector<pair<double, vector<vec3_t>>> read_eigs();



void plot_modes(string filename, airebo::system_control_t &sys,
    vector<pair<double, vector<vec3_t>>> modes);

void write_spectra(vector<pair<double, vector<vec3_t>>> modes,
    structure_t structure, const int *pol_axes, const string &filename);

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
        if (generate_eigs(settings.phonopy_settings) != 0)
            return EXIT_FAILURE;

        // read eigenmodes/etc
        auto modes = read_eigs();

        cout << "num modes: " << modes.size() << endl;
        if (!modes.empty())
            cout << "num eigs: " << modes[0].second.size() << endl;

        airebo::system_control_t sys;
        sys.init(structure);
        plot_modes("modes.dat", sys, modes);

        write_spectra(modes, structure,
            settings.phonopy_settings.polarization_axes, "spectra.dat");

        int xx[2] = {0, 0},
            yy[2] = {1, 1},
            zz[2] = {2, 2},
            xy[2] = {0, 1};
        write_spectra(modes, structure, xx, "spectra_xx.dat");
        write_spectra(modes, structure, yy, "spectra_yy.dat");
        write_spectra(modes, structure, zz, "spectra_zz.dat");
        write_spectra(modes, structure, xy, "spectra_xy.dat");

        write_spectra(modes, structure, nullptr, "spectra_avg.dat");
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
    vector<double> forces = gradient;
    vscal(-1.0, forces);

    outfile << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"
            << "<varray name=\"forces\" >\n";

    for (int i = 0; i < forces.size(); i += 3)
    {
        outfile << "<v> ";
        for (int j = 0; j < 3; ++j)
            outfile << forces[i + j] << " ";

        outfile << "</v>\n";
    }

    outfile << "</varray>" << endl;
}

void relax_structure(structure_t &structure, run_settings_t rset)
{
    auto supercell = util::construct_supercell(structure,
        rset.phonopy_settings.supercell_dim);

    // construct system + minimize with default parameters
    if (rset.add_hydrogen)
    {
        airebo::system_control_t sys;
        sys.init(supercell);
        sys.add_hydrogen();
        supercell = sys.get_structure();
    }

//    airebo::system_control_t sys;
    lammps::system_control_t sys;
    sys.init(supercell, rset.lammps_settings);

    minimize::acgsd_settings_t settings;
    settings.gradient_tolerance = 1e-5;

    minimize::acgsd(sys.get_diff_fn(), sys.get_position(), settings);

    supercell = sys.get_structure();

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
    std::unique_ptr<lammps::system_control_t> sys;
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

//        airebo::system_control_t sys;
        if (i == 1)
        {
            sys = std::make_unique<lammps::system_control_t>();
            sys->init(structure, rset.lammps_settings);
        }
        else
        {
            sys->set_position(structure.positions);
        }

        sys->update();

        output_xml(output_file, sys->get_gradient());

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


void write_spectra(vector<pair<double, vector<vec3_t>>> modes,
    structure_t structure, const int *pol_axes, const string &filename)
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
    if (pol_axes != nullptr)
    {
        spectra = sp2::phonopy::raman_spectra(pol_vecs[pol_axes[0]],
            pol_vecs[pol_axes[1]], temperature, modes, structure);
    }
    else
    {
        spectra = sp2::phonopy::raman_spectra_avg(true, temperature,
            modes, structure);
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

    ofstream outfile(filename);
    for (auto &shift : spectra)
    {
        outfile << shift.first << ' '
                << shift.second / maxi << ' '
                << shift.second << endl;

        for (std::size_t i = 0; i < bins.size(); ++i)
            bins[i] += gaussian(i * bin_step, shift.first) * shift.second;
    }
    outfile.close();

    outfile.open("gauss_" + filename);
    double max_bin = *std::max_element(bins.begin(), bins.end());
    for (std::size_t i = 0; i < bins.size(); ++i)
        outfile << (i * bin_step) << ' '
                << bins[i] / max_bin << std::endl;

    outfile.close();
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
