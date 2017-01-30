#include "common/io/util.hpp"
#include "common/io/file_types/xyz.hpp"

#include <sstream>
#include <iostream>
#include <fstream>
#include <map>

using namespace std;
using namespace sp2;

/// move the input stream to the next xyz frame
bool next_frame(ifstream& infile)
{
    string line;
    if (!getline(infile, line) || line.empty())
        return false;

    // read in the number of atoms in the frame
    int num_atoms;
    try {
        num_atoms = stoi(line);
    } catch (std::invalid_argument &e) {
        return false;
    }

    // just discard all the lines until the next frame
    while (num_atoms >= 0 && getline(infile, line))
        num_atoms -= 1;

    // check if all lines were discarded properly,
    // num_atoms should be -1 (not zero) because of
    // the extra comment line
    return (num_atoms == -1);
}

/// read a frame of an xyz file
bool io::read_xyz(std::string filename, std::string &comment_line,
    std::vector<std::string> &types, std::vector<double> &positions,
    int frame)
{
    const string info = "io::read_xyz(): filename: \"" + filename + "\"";

    ifstream infile(filename);
    if (!infile.is_open())
    {
        cout << "Error opening file. " << info << endl;
        return false;
    }

    int current_frame = 0;
    while (current_frame < frame && next_frame(infile))
        ++current_frame;

    if (current_frame != frame)
    {
        cout << "Failed to read frame " << frame << ". "
             << info << endl;
        return false;
    }

    // read in the number of atoms in the frame
    string line;
    if (!getline(infile, line))
        return false;

    int num_atoms;
    try {
        num_atoms = stoi(line);
    } catch (std::invalid_argument &e) {
        cout << "Failed to parse number of atoms \"" <<
             line << "\". " << info << endl;
        return false;
    }

    comment_line.clear();
    if (!getline(infile, comment_line))
        return false;

    types.clear();
    positions.clear();

    while (num_atoms > 0 && getline(infile, line))
    {
        num_atoms -= 1;
        if (line.empty())
            return false;

        vector<string> split = io::split_line(line, " \t");
        if (split.size() < 4)
            return false;

        types.push_back(split[0]);
        for (int i = 0; i < 3; ++i)
            positions.push_back(stod(split[i + 1]));
    }

    if (num_atoms != 0)
    {
        cout << "Error, read less than specified number of atoms. "
             << info << endl;
        return false;
    }

    return true;
}

bool io::write_xyz(std::string filename, const std::string &comment_line,
    const std::vector<std::string> &types, const std::vector<double> &positions,
    bool append)
{
    const string info = "io::write_xyz(): filename: \"" + filename + "\"";

    if (types.size() != positions.size() / 3)
    {
        cout << "Error, mismatch between atom type and position vector size. "
             << info << endl;
        return false;
    }

    ofstream outfile;
    if (append)
        outfile.open(filename, fstream::app);
    else
        outfile.open(filename);

    if (!outfile.is_open())
    {
        cout << "Error opening file. " << info << endl;
        return false;
    }

    outfile.precision(8);

    std::size_t num_atoms = types.size();
    outfile << num_atoms << '\n'
            << comment_line << '\n';

    for (std::size_t i = 0; i < num_atoms; ++i)
        outfile << types[i] << '\t'
                << positions[i * 3] << '\t'
                << positions[i * 3 + 1] << '\t'
                << positions[i * 3 + 2] << '\n';

    return true;
}

/// get the number of frames in an xyz file
int get_xyz_frames(std::string filename)
{
    ifstream infile(filename);
    if (!infile.is_open())
        return -1;

    int frames = 0;
    while (next_frame(infile))
        ++frames;

    return frames;
}
