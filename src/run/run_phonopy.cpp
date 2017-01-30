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


using namespace std;
using namespace sp2;


int csyscall_fork(string command, char *args[])
{
    int pid = fork();

    if (pid == 0)
    {
        prctl(PR_SET_PDEATHSIG, SIGKILL);
        execvp(command.c_str(), args);

        return EXIT_FAILURE;
    }
    else if (pid > 0)
    {
        int status = EXIT_FAILURE;
        waitpid(pid, &status, 0);

        return status;
    }
    else
    {
        return pid;
    }
}

int syscall_fork(string command, vector<string> args)
{
    vector<vector<char>> cargs;

    args.insert(args.begin(), command);

    vector<char*> cargs_actual;
    for (auto s : args)
    {
        cargs.emplace_back();
        vector<char> &temp = cargs.back();
        for (auto c : s)
            temp.push_back(c);

        temp.push_back('\0');

        cargs_actual.push_back(&temp[0]);
    }

    cargs_actual.push_back(nullptr);

    return csyscall_fork(command, &cargs_actual[0]);
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

/// read all input poscar files, return vector of output filenames
vector<string> process_displacements()
{
    auto get_ext = [](int i) {
        stringstream ss;
        ss << setw(3) << setfill('0') << i;
        return ss.str();
    };

    vector<string> output_filenames;
    for (int i = 1; i < 1000; ++i)
    {
        string input_file = "POSCAR-" + get_ext(i),
            output_file = "vasprun.xml-" + get_ext(i);

        structure_t structure;
        if (!io::file_exists(input_file) ||
            !io::read_structure(input_file, structure, file_type::POSCAR))
        {
            break;
        }

        airebo::system_control_t sys;
        sys.init(structure);

        sys.update();
        output_xml(output_file, sys.get_gradient());

        output_filenames.push_back(output_file);
    }

    return output_filenames;
}

constexpr double kronecker_delta(size_t a, size_t b)
{
    return a == b ? 1 : 0;
}

double polarization(vector<vec3_t> &eigs, airebo::system_control_t &sys,
    int A, int B)
{
    constexpr double const1 = 0.32, //  a|| -   a|- in Angstroms^3
        const2 = 7.55,              // a'|| + 2a'|- in Angstroms^2
        const3 = 2.60;              // a'|| -  a'|- in Angstroms^2

#warning error with get_graph() for repeated bonds across periodic bound
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

#warning may need to scale with mass
        vec3_t eig = eigs[bond.a],// * (types[bond.a] == atom_type::HYDROGEN ? 1 : 1 / sqrt(12)),
            unit_delta = bonds[bond.id].unit_vector();

        double term1 = (const2 / 3) * dot(unit_delta, eig) * kronecker_delta(A, B),
            term2 = const3 * (unit_delta[A] * unit_delta[B] - kronecker_delta(A, B) / 3) * dot(unit_delta, eig),
#warning might be R[A] X[B] PLUS R[B] X[A]
            term3 = (const1 / len) * (unit_delta[A] * eig[B] + unit_delta[B] * eig[A] - 2 * unit_delta[A] * unit_delta[B] * dot(unit_delta, eig));

        double alternate = term1 + term2 + term3;

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

#warning assumes polarization unit vectors are on an axis
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

void relax_structure(structure_t &structure, bool add_hydrogen)
{
    // construct system + minimize with default parameters
    airebo::system_control_t sys;
    sys.init(structure);
    if (add_hydrogen)
        sys.add_hydrogen();

    minimize::acgsd_settings_t settings;
    settings.gradient_tolerance = 1e-5;

    minimize::acgsd(sys.get_diff_fn(), sys.get_position(), settings);

//    double min_ptnl = sys.get_value();
//    auto min_structure = sys.get_structure();
//    double min_x = 0;
//
//    for (double x = 0; x < 0.01; x += 0.000001)
//    {
//        double lattice[3][3] = {
//            {4.254 - x, 0, 0},
//            {0, 25, 0},
//            {0, 0, 25}
//        };
//        sys.set_lattice(lattice);
//        sys.update();
//
//        minimize::acgsd(sys.get_diff_fn(), sys.get_position(), settings);
//        if (sys.get_value() < min_ptnl)
//        {
//            min_ptnl = sys.get_value();
//            min_structure = sys.get_structure();
//            min_x = x;
//        }
//    }
//
//    cout << "min ptnl: " << min_ptnl << ", min lattice: " << min_structure.lattice[0][0]
//         << ", min x: " << min_x << endl;

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

    return syscall_fork("phonopy", {"-d", "disp.conf"});
}

int generate_force_sets()
{
    // process displacements
    auto force_filenames = process_displacements();

    // generate force sets, "phonopy -f file1 file2 ..."
    vector<string> args = {"-f"};
    for (auto file : force_filenames)
        args.push_back(file);

    return syscall_fork("phonopy", args);
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
            << "FORCE_CONSTANTS = WRITE\n";

    outfile.close();

    outfile.open("band_labels.txt");
    for (auto qp : pset.qpoints)
        outfile << '\t' << qp.label;

    outfile.close();

    return syscall_fork("phonopy", {"band.conf"});
}

int generate_eigs(phonopy::phonopy_settings_t pset)
{
    ofstream outfile("eigs.conf");

    outfile << "DIM = " << pset.supercell_dim[0]
                << ' ' << pset.supercell_dim[1]
                << ' ' << pset.supercell_dim[2] << '\n'
            << "BAND = 0 0 0   1/2 0 0\n"
            << "BAND_POINTS = 1\n"
            << "FORCE_CONSTANTS = WRITE\n"
            << "EIGENVECTORS = .TRUE.\n";

    outfile.close();

    return syscall_fork("phonopy", {"eigs.conf"});
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

    return syscall_fork("phonopy", {"dos.conf"});
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

void remove_overlap(structure_t &structure)
{
    auto positions = dtov3(structure.positions);
    cout << positions.size() << endl;

    vec3_t lattice_vec[3] = {
        vec3_t(structure.lattice[0]),
        vec3_t(structure.lattice[1]),
        vec3_t(structure.lattice[2])
    };

    auto get_offset = [&lattice_vec](int i, int j, int k) {
        return i * lattice_vec[0] +
               j * lattice_vec[1] +
               k * lattice_vec[2];
    };

    auto gen_positions = [&](){
        vector<vec3_t> pos_all;
        for (auto v : positions)
            for (int i = -1; i <= 1; ++i)
                for (int j = -1; j <= 1; ++j)
                    for (int k = -1; k <= 1; ++k)
                        if (i != 0 || j != 0 || k != 0)
                            pos_all.push_back(v + get_offset(i, j, k));

        return pos_all;
    };

    do {
        auto pos_all = gen_positions();

        bool ok = true;
        for (auto it = positions.begin(); it != positions.end(); ++it)
        {
            auto pos_a = *it;
            for (auto pos_b : pos_all)
            {
                if ((pos_b - pos_a).mag() < 1e-3)
                {
                    ok = false;
                    break;
                }
            }

            if (!ok)
            {
                positions.erase(it);
                break;
            }
        }

        if (ok)
            break;

    } while (true);

    structure.positions = v3tod(positions);
    structure.types.resize(positions.size());
    cout << positions.size() << endl;
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

void find_breathing_mode(structure_t structure,
    vector<pair<double, vector<vec3_t>>> &modes, int last_mode)
{
    constexpr double tol = 0.01;
    constexpr int mode_tol = 5;
    const int mode_target = last_mode + 10;

    double middle = 3.684072;
    auto pos = dtov3(structure.positions);

    vec3_t avg;
    for (auto v : pos)
        avg += v / pos.size();

    auto fitness_fn = [=](double dy, double ey) {
        auto m = dy * ey;
        return m > 0 ? m : 0;
    };

    double max_sum = -1e9;
    int id = 0, max_id = 0;
    pair<double, vector<vec3_t>> breathing_mode;

    for (int id = mode_target - mode_tol; id < mode_target + mode_tol; ++id)
    {
        auto &eigs = modes[id].second;

        vec3_t low{},
            high{};
        for (int i = 0; i < pos.size(); ++i)
        {
            if (structure.types[i] == atom_type::HYDROGEN)
                continue;

            if (pos[i].y() > avg.y())
                high += eigs[i];
            else
                low += eigs[i];
        }

        double val = dot(high, -low);
        if (val > max_sum)
        {
            breathing_mode = modes[id];
            max_sum = val;
            max_id = id;
        }
    }

    ofstream outfile("bm.dat");
    outfile << max_id << ' ' << breathing_mode.first << ' ' << max_sum << endl;
    outfile.close();

    cout << max_id << ' ' << breathing_mode.first << ' ' << max_sum << endl;
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
            {
                outfile << pos1[k] << ' ';
            }

            for (int k = 0; k < 3; ++k)
            {
                outfile << pos2[k] << ' ';
            }

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
    stringstream ss;
    ss << "{\"status\":\"" << desc << "\"}";

    string string = ss.str();
    ofstream outfile(filename);
    outfile.write(string.c_str(), string.size());
    outfile.close();
}

int sp2::run_phonopy(const run_settings_t &settings, MPI_Comm comm)
{
    auto structure = settings.structure;

    // relax
    write_log(settings.log_filename, "Relaxing Structure");
    relax_structure(structure, settings.add_hydrogen);

    sort_structure_types(structure);
    io::write_structure("relaxed.xyz", structure);

    // write poscar
    io::write_structure("POSCAR", structure, false, file_type::POSCAR);

    write_log(settings.log_filename, "Generating Displacements");
    if (generate_displacements(settings.phonopy_settings) != 0)
        return EXIT_FAILURE;

    write_log(settings.log_filename, "Generating Force Sets");

    if (generate_force_sets() != 0)
        return EXIT_FAILURE;

//    generate_eigs(dim)

    write_log(settings.log_filename, "Computing Force Constants");
    if (generate_bands(settings.phonopy_settings) != 0)
        return EXIT_FAILURE;

//    // read eigenmodes/etc
//    auto modes = read_eigs();
////    find_breathing_mode(structure, modes, last_mode);
//
//    cout << "num modes: " << modes.size() << endl;
//    if (!modes.empty())
//        cout << "num eigs: " << modes[0].second.size() << endl;
//
//    airebo::system_control_t sys;
//    sys.init(structure);
//    plot_modes("modes.dat", sys, modes);
//
//    int pol_axes[2] = {0, 0};
//    io::deserialize_array(config["polarization_axes"], &pol_axes[0]);
//    cout << "pol_axes " << pol_axes[0] << ' ' << pol_axes[1] << endl;
//
//    write_spectra(modes, sys, pol_axes, "spectra.dat");
//
//
//    int xx[2] = {0, 0},
//        yy[2] = {1, 1},
//        zz[2] = {2, 2},
//        xy[2] = {0, 1};
//    write_spectra(modes, sys, xx, "spectra_xx.dat");
//    write_spectra(modes, sys, yy, "spectra_yy.dat");
//    write_spectra(modes, sys, zz, "spectra_zz.dat");
//    write_spectra(modes, sys, xy, "spectra_xy.dat");

    return EXIT_SUCCESS;
}
