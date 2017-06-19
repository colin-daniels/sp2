#include "common/io/util.hpp"

#include <iostream>
#include <sstream>
#include <fstream>

using namespace std;
using namespace sp2;

bool io::clear_file(std::string filename)
{
    ofstream outfile(filename.c_str());
    if (!outfile.is_open())
    {
        cerr << "Error opening " << filename << " in clear_file()!" << endl;
        return false;
    }
    outfile.close();
    return true;
}

bool io::remove_file(std::string filename)
{
    if (filename.empty() || !file_exists(filename))
        return false;
    return remove(filename.c_str()) == 0;
}

bool io::file_exists(std::string filename)
{
    ifstream infile(filename.c_str());
    if (!infile.is_open())
        return false;

    infile.close();
    return true;
}

bool io::read_file(std::string filename, std::string &content)
{
    ifstream infile(filename, fstream::binary);
    if (!infile.is_open())
        return false;

    infile.seekg(0, fstream::end);
    content.resize(infile.tellg(), '\0');

    infile.seekg(0, fstream::beg);
    infile.read(&content[0], content.size());

    return true;
}

bool io::write_file(std::string filename, const std::string &content)
{
    ofstream outfile(filename, fstream::binary);
    if (!outfile.is_open())
        return false;

    outfile.write(&content[0], content.size());
    return true;
}

bool io::copy_file(std::string source, std::string dest)
{
    if (dest.empty() || source.empty())
    {
        cerr << "Error, source or destination string is empty in copy_file" << endl;
        return false;
    }
    if (source[source.size() - 1] == '/' || source[source.size() - 1] == '\\')
    {
        cerr << "Error, source for copy is a directory: " << source << endl;
        return false;
    }

    if (dest[dest.size() - 1] == '/' || dest[dest.size() - 1] == '\\')
    {
        size_t pos_back = source.find_last_of("/\\");
        if (pos_back == string::npos)
            dest = dest + source;
        else
            dest = dest + source.substr(pos_back + 1);
    }

    ifstream infile(source.c_str(), fstream::binary);
    ofstream outfile(dest.c_str(), fstream::binary);

    if (!infile.is_open())
    {
        cerr << "failed opening input file for copy " << source << endl;
        return false;
    }
    else if (!outfile.is_open())
    {
        cerr << "failed opening output file for copy " << dest << endl;
        return false;
    }

    outfile << infile.rdbuf();

    outfile.close();
    infile.close();
    return true;
}


void io::trim(std::string &str, std::string chars)
{
    size_t left = str.find_first_not_of(chars, 0);
    size_t right = str.find_last_not_of(chars, str.size());

    // if no non-whitespace character was found, delete the string content and return
    if (left == string::npos)
    {
        str.clear();
        return;
    }
    str = str.substr(left, right - left + 1);
}

std::vector<std::string> io::split_line(const std::string &line,
    const std::string &chars)
{
    vector<string> output;
    string temp;

    for (unsigned int i = 0; i < line.size(); ++i)
    {
        char c = line[i];

        // if the current character was specified to be split on
        if (chars.find(c) != string::npos)
        {
            if (temp.empty())
                continue;
            else
            {
                output.push_back(temp);
                temp.clear();
            }
        }
        else // if it is not
            temp.push_back(c);
    }
    if (!temp.empty())
        output.push_back(temp);
    return output;
}

std::string io::toupper(const std::string &input)
{
    std::string output = input;
    for (auto &c : output)
        c = static_cast<char>(::toupper(c));
    return output;
}

std::string io::tolower(const std::string &input)
{
    std::string output = input;
    for (char &c : output)
        c = static_cast<char>(::tolower(c));
    return output;
}

bool sp2::io::move_file(std::string source, std::string dest)
{
    return io::copy_file(source, dest) && io::remove_file(source);
}
