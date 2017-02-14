#include "run/run_types.hpp"
#include "common/io/util.hpp"
#include "common/io/structure.hpp"
#include "airebo/system_control_t.hpp"
#include "common/util/blas.hpp"
#include "common/minimize/minimize.hpp"
#include "common/vec3_t.hpp"
#include "common/json/json.hpp"

#include <cstdlib>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <common/vec3_t.hpp>
#include <iostream>
#include <cstdlib>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/prctl.h>

#ifdef SP2_ENABLE_LAMMPS
#include "lammps/lammps_interface.hpp"
#endif

using namespace std;
using namespace sp2;

double polarization(vector<vec3_t> &eigs, airebo::system_control_t &sys,
    int A, int B);

vector<pair<double, double>> calc_raman_spectra(
    vector<pair<double, vector<vec3_t>>> modes,
    airebo::system_control_t &sys, double temperature, double incident,
    vec3_t pol_A, vec3_t pol_B);

void relax_structure(structure_t &structure, run_settings_t rset);

int generate_displacements(phonopy::phonopy_settings_t pset);
int generate_force_sets(run_settings_t rset);

int generate_bands(phonopy::phonopy_settings_t pset);
int generate_eigs(phonopy::phonopy_settings_t pset);
int generate_dos(phonopy::phonopy_settings_t pset);

vector<pair<double, vector<vec3_t>>> read_eigs();
void repeat_structure(structure_t &structure, int n_times);


void plot_modes(string filename, airebo::system_control_t &sys,
    vector<pair<double, vector<vec3_t>>> modes);

void write_spectra(const vector<pair<double, vector<vec3_t>>> &modes,
    airebo::system_control_t &sys, const int *pol_axes, const string &filename);

void write_log(string filename, string desc);

void output_xml(string filename, const vector<double> &gradient);


int sp2::run_phonopy(const run_settings_t &settings, MPI_Comm)
{
    auto structure = settings.structure;
    io::write_structure("input.xyz", structure);

    // relax
    write_log(settings.log_filename, "Relaxing Structure");
    relax_structure(structure, settings);

    sort_structure_types(structure);
    io::write_structure("relaxed.xyz", structure);

    io::write_structure("POSCAR", structure, false, file_type::POSCAR);

    write_log(settings.log_filename, "Generating Displacements");
    if (generate_displacements(settings.phonopy_settings) != 0)
        return EXIT_FAILURE;

    write_log(settings.log_filename, "Generating Force Sets");

    if (generate_force_sets(settings) != 0)
        return EXIT_FAILURE;

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

        write_spectra(modes, sys, settings.phonopy_settings.polarization_axes,
            "spectra.dat");
    }

    if (settings.phonopy_settings.calc_bands &&
        generate_bands(settings.phonopy_settings) != 0)
        return EXIT_FAILURE;

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

constexpr double kronecker_delta(int a, int b)
{
    return a == b ? 1 : 0;
}

double polarization(vector<vec3_t> &eigs, airebo::system_control_t &sys,
    int A, int B)
{
    constexpr double const1 = 0.32, //  a|| -   a|- in Angstroms^3
        const2 = 7.55,              // a'|| + 2a'|- in Angstroms^2
        const3 = 2.60;              // a'|| -  a'|- in Angstroms^2

// TODO: error with get_graph() for repeated bonds across periodic bound
    auto graph = sys.get_bond_control().get_graph();
    auto bonds = dtov3(sys.get_bond_control().get_bond_deltas());
    auto types = sys.get_types();

    double test_sum = 0;
    for (auto e1 : eigs)
        test_sum += dot(e1, e1);

//    cout << "eigs dot: " << test_sum << endl;

    double sum = 0;
    for (graph::ud_edge_t bond : graph.edges())
    {
        double len = bonds[bond.id].mag();
        if (types[bond.a] == atom_type::HYDROGEN ||
            types[bond.b] == atom_type::HYDROGEN)
            continue;

// TODO: may need to scale with mass
        vec3_t eig = eigs[bond.a],// * (types[bond.a] == atom_type::HYDROGEN ? 1 : 1 / sqrt(12)),
            unit_delta = bonds[bond.id].unit_vector();

        double term1 = (const2 / 3) * dot(unit_delta, eig)
                       * kronecker_delta(A, B),
            term2 = const3 * (unit_delta[A] * unit_delta[B]
                      - kronecker_delta(A, B) / 3) * dot(unit_delta, eig),
// TODO: might be R[A] X[B] PLUS R[B] X[A]
            term3 = (const1 / len)  * (unit_delta[A] * eig[B]
                    + unit_delta[B] * eig[A]
                    - 2 * unit_delta[A] * unit_delta[B] * dot(unit_delta, eig));

        sum += term1 + term2 + term3;
    }

    return -sum;
}

vector<pair<double, double>> calc_raman_spectra(
    vector<pair<double, vector<vec3_t>>> modes,
    airebo::system_control_t &sys, double temperature, double incident,
    vec3_t pol_A, vec3_t pol_B)
{
    vector<pair<double, double>> result;

    for (auto mode : modes)
    {
        double frequency = mode.first;
        vector<vec3_t> &eigs = mode.second;

        // hbar / k_b = 7.63823×10^-12 s K
        const double hkt = 7.63823e-12 / temperature;
        // thermal average occupation number for this mode/temperature
        // <n(omega_f)> = [exp(hbar omega_f / k_b T) - 1]^-1
        double n_omega = 1.0 / (exp(hkt * frequency) - 1.0);
        double prefactor = (n_omega + 1) / frequency;

// TODO: assumes polarization unit vectors are on an axis
        double pol = 0;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                pol += pol_A[i] * pol_B[j] *
                    polarization(eigs, sys, i, j);

        const double scattered = incident - frequency;
        double intensity = /*incident * (scattered * scattered * scattered)
             * prefactor  */ (pol * pol);

        // frequency is in Hz so we need to convert to cm^-1
        double wavenumber = frequency * 33.35641e-12;
        result.emplace_back(wavenumber, intensity);
    }

    return result;
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

//    airebo::system_control_t sys;
    lammps::system_control_t sys;
    sys.init(structure, rset.lammps_settings);

    minimize::acgsd_settings_t settings;
    settings.gradient_tolerance = 1e-5;

    minimize::acgsd(sys.get_diff_fn(), sys.get_position(), settings);

    structure = sys.get_structure();
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

            // phonopy outputs frequency in terahertz
            frequency *= 1e12;

            modes.emplace_back(frequency, vector<vec3_t>{});
            reading_mode = true;
        }
    }

    return modes;
}

void repeat_structure(structure_t &structure, int n_times)
{
    vector<double> new_pos;
    vector<atom_type> new_types;

    for (int i = 0; i < n_times; ++i)
    {
        for (int j = 0; j < structure.positions.size(); ++j)
        {
            new_pos.push_back(structure.positions[j]);
            if (j % 3 == 0)
                new_pos.back() += i * structure.lattice[0][0];

        }
        copy(structure.types.begin(), structure.types.end(),
            back_inserter(new_types));
    }

    structure.lattice[0][0] *= n_times;
    structure.positions = new_pos;
    structure.types = new_types;
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


void write_spectra(const vector<pair<double, vector<vec3_t>>> &modes,
    airebo::system_control_t &sys, const int *pol_axes, const string &filename)
{
    double incident_freq = 2.818e14, // 1064 nm
        temperature = 293;

    vec3_t pol_vecs[3] = {
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1}
    };

    auto spectra = calc_raman_spectra(modes, sys, temperature, incident_freq,
        pol_vecs[pol_axes[0]], pol_vecs[pol_axes[1]]);

    double maxi = 0;
    for (auto mode : spectra)
        maxi = max(mode.second, maxi);

    ofstream outfile(filename);
    for (auto &shift : spectra)
        outfile << shift.first << ' '
                << shift.second / maxi << ' '
                << shift.second << endl;
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
