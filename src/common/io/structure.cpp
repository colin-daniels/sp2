#include "common/io/structure.hpp"
#include "common/io/util.hpp"

#include "common/json/json.hpp"
#include "common/io/file_types/xyz.hpp"
#include "common/io/file_types/cif.hpp"
#include "common/io/file_types/poscar.hpp"

#include <iostream>
#include <unordered_map>

using namespace std;
using namespace sp2;

file_type io::infer_filetype(std::string filename)
{
    static const unordered_map<string, file_type> type_map = {
        {"xyz",    file_type::XYZ},
        {"json",   file_type::JSON},
        {"cif",    file_type::CIF},
        {"cifs",    file_type::CIF},
        {"vasp",   file_type::POSCAR},
        {"poscar",   file_type::POSCAR}
    };

    // if it looks like a filename, pull the extension off
    string::size_type pos = filename.find_last_of('.');
    if (pos != string::npos && pos + 1 < filename.size())
        filename = filename.substr(pos + 1);

    for (auto &c : filename)
        c = std::tolower(c);

    return type_map.count(filename) ?
           type_map.at(filename) : file_type::UNKNOWN;
}

/// write structure to a file
bool io::write_structure(std::string filename,
        const structure_t& input, bool append, file_type type)
{
    if (type == file_type::AUTO || type == file_type::UNKNOWN)
        type = infer_filetype(filename);

    switch (type)
    {
    case file_type::XYZ: {
        Json::Value val;
        val["lattice"] = input.serialize_lattice();

        // convert types to strings
        vector<std::string> types;
        for (atom_type t : input.types)
            types.push_back(enum_to_str<atom_type>(t));

        // write the lattice (unstyled, no line breaks) to the comment line
        // and then write the file
        string comment_line;
        return io::write_json(val, comment_line, false) &&
               io::write_xyz(filename, comment_line,
                       types, sp2::v3tod(input.positions), append);
    }
    case file_type::JSON: {
        // note: no append option at the moment

        // just serialize to json and write
        Json::Value val;
        return input.serialize(val) &&
               io::write_json_file(val, filename);
    }
    case file_type::CIF:
        return io::write_cif(filename, input);
    case file_type::POSCAR:
        return io::write_poscar(filename, input);
    default:
        cerr << "Failed to write structure file \"" + filename + "\""
             << " with unknown file type.";
        return false;
    }
}

/// read structure from a file
bool io::read_structure(std::string filename,
        structure_t &output, file_type type)
{
    if (type == file_type::AUTO || type == file_type::UNKNOWN)
        type = infer_filetype(filename);

    switch (type)
    {
    case file_type::XYZ: {
        // parse the types and positions normally
        string comment_line;
        vector<std::string> types;

        std::vector<double> pos_in;
        if (!io::read_xyz(filename, comment_line, types, pos_in))
            return false;

        output.positions = sp2::dtov3(pos_in);

        // convert strings to types
        for (string t : types)
            output.types.push_back(enum_from_str<atom_type>(t));

        // and then deserialize the lattice information from
        // the comment line
        Json::Value val;
        if (Json::Reader().parse(comment_line, val))
            return output.deserialize_lattice(val["lattice"]);
        else
        {
            if (!comment_line.empty())
                std::cerr << "Warning, non-empty comment failed "
                          << "to parse into JSON." << std::endl;

            return true;
        }
    }
    case file_type::JSON: {
        // just read and deserialize
        Json::Value val;
        return io::read_json_file(val, filename) &&
               output.deserialize(val);
    }
    case file_type::POSCAR:
        return io::read_poscar(filename, output);
    default:
        cerr << "Failed to read structure file \"" + filename + "\""
             << " with unknown file type.";
        return false;
    }
}

bool io::convert_file(std::string input_file, std::string output_file,
        file_type input_type, file_type output_type)
{
    // interchange uses the structure_t type
    structure_t structure;
    return read_structure(input_file, structure, input_type) &&
           write_structure(output_file, structure, false, output_type);
}
