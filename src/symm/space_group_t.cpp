#include "symm/space_group_t.hpp"
#include "common/io/util.hpp"
#include "common/util/templates.hpp"

#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <regex>
#include <cmath>

using namespace std;
using namespace sp2;

/// parse a single line into transformation matrix
bool parse_transform(string line, array<double, 16> &output);

int symm::space_group_t::get_number() const {
    return number;}

std::size_t symm::space_group_t::n_symm() const {
    return transform.size();}

std::string symm::space_group_t::get_name() const {
    return name;}

std::vector<vec3_t> symm::space_group_t::apply_symm(
    const std::vector<vec3_t> &input) const
{
    // make sure the output can hold everything
    std::vector<vec3_t> output;
    output.reserve(input.size() * transform.size());

    // apply transformations to every input vector in order
    for (auto &mat : transform)
    {
        for (auto &vec : input)
        {
            auto new_vec = vec3_t();
            for (auto i = 0; i < 3; ++i)
                new_vec[i] = mat[i * 4] * vec.x()
                             + mat[i * 4 + 1] * vec.y()
                             + mat[i * 4 + 2] * vec.z()
                             + mat[i * 4 + 3];

            output.push_back(new_vec);
        }
    }

    return output;
}

std::istream& symm::space_group_t::parse(std::istream &stream)
{
    if (!stream)
        return stream;

    string line;
    transform.clear();

    for (int line_n = 0; getline(stream, line) && !line.empty(); ++line_n)
    {
        if (line_n == 0)
            number = stoi(line);
        else if (line_n == 1)
            line.clear(); // ignore it
        else if (line_n == 2)
            name = line;
        else
        {
            transform.push_back(array<double, 16>());
            if (!parse_transform(line, transform.back()))
                throw runtime_error("invalid space group input file, line "
                                    + to_string(line_n) + ": '" + line + "'.");

            // if (name == "I a -3 d")
            // {
            //     cout << line << endl;
            //     for (size_t i = 0; i < 4; ++i)
            //     {
            //         for (size_t j = 0; j < 4; ++j)
            //             cout << transform.back()[i * 4 + j] << ' ';
            //         cout << endl;
            //     }
            // }
        }
        line.clear();
    }

    return stream;
}

bool parse_transform(string line, array<double, 16> &output)
{
    // cout << line << endl;
    static const regex offset_regex("([+-]?\\d)/(\\d)"),
        coordinate_regex("[+-]?[xyz]");

    // make sure output is cleared
    for (auto &d : output) {d = 0;}

    // note: row looks like "1/2+y,1/2+x,-z"
    // so first split on commas
    auto columns = io::split_line(line, ",");

    // check size
    if (columns.size() != 3)
        return false;

    // helper lambda for getting matches in the form of a range
    auto matches = [](const std::string &str, const std::regex &reg) {
        return make_range(
            sregex_iterator(str.begin(), str.end(), reg),
            sregex_iterator()
        );
    };

    for (std::size_t i = 0; i < 3; ++i)
    {
        auto column = columns[i];

        // find any offsets
        for (auto offset : matches(column, offset_regex))
            output[i * 4 + 3] += stod(offset[1].str())
                                 / stod(offset[2].str());

        // find coordinates
        for (auto coord : matches(column, coordinate_regex))
        {
            auto str  = coord.str();
            // x,y,z -> 0,1,2
            auto axis = str.back() - 'x';

            output[i * 4 + axis] += (str.front() == '-') ? -1 : 1;
        }
    }

    // last row of the matrix is just (0, 0, 0, 1)
    output[15] = 1;

    // // make sure the line is valid
    // if (distance(matches.begin(), matches.end()) != 3)
    //     return false;

    // // fill each row of the 4x4 transformation matrix
    // size_t row = 0;
    // for (auto match : matches)
    // {
    //     // if there is a translation (e.g. 1/4)
    //     auto trans_str = match[1].str();
    //     if (!trans_str.empty())
    //         output[row * 4 + 3] = stod(trans_str.substr(0, 1))
    //                             / stod(trans_str.substr(2, 1));

    //     // parse axis (or axes in the case of something like 1/2+x-y)
    //     for (int i : {2, 3})
    //     {
    //         auto axis_str = match[i].str();
    //         if (axis_str.empty())
    //             continue;

    //         auto axis = axis_str.back() - 'x';
    //         output[row * 4 + axis] = 1;

    //         if (axis_str[0] == '-')
    //             output[row * 4 + axis] = -1;
    //     }

    //     // increment the row
    //     row += 1;
    // }

    // // last row of the matrix is just (0, 0, 0, 1)
    // output[row * 4 + 3] = 1;

    return true;
}

std::map<std::string, symm::space_group_t>
symm::read_groups(std::string filename)
{
    map<string, space_group_t> groups;

    ifstream infile(filename);
    if (!infile.is_open())
    {
        cerr << "Failed to open space group input file \""
             << filename << "\"." << endl;
        return groups;
    }

    while (infile)
    {
        space_group_t group;
        group.parse(infile);

        string name = group.get_name();
        if (!name.empty())
            groups[name] = group;
    }

    return groups;
}
