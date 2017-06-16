#include "common/json/json.hpp"
#include "common/io/util.hpp"

#include <iostream>
#include <string>

using namespace std;
using namespace sp2;

bool io::write_json(const Json::Value &root,
    std::string &output, bool styled)
{
    if (!root)
        return false;

    if (styled)
        output = Json::StyledWriter().write(root);
    else
    {
        Json::FastWriter writer;
        writer.omitEndingLineFeed();
        output = writer.write(root);
    }
    return true;
}

bool io::write_json_file(const Json::Value &root,
    std::string filename, bool styled)
{
    string content;
    if (!io::write_json(root, content, styled))
        return false;

    // add an extra line feed if set to not styled, since this is an output file
    if (!styled)
        content += "\n";

    if (!io::write_file(filename, content))
    {
        cerr << "Failed to write JSON file \"" << filename << "\"" << endl;
        return false;
    }
    return true;
}

bool io::read_json(Json::Value &root, const std::string &input)
{
    Json::Reader reader;
    if (!reader.parse(input, root))
    {
        cerr << "Failed to parse JSON\n"
             << reader.getFormattedErrorMessages();
        return false;
    }

    return true;
}

bool io::read_json_file(Json::Value &root, std::string filename)
{
    // load the file into a string
    string content;
    if (!io::read_file(filename, content))
    {
        cerr << "Failed to open (or empty) JSON file \""
             << filename << "\"" << endl;
        return false;
    }

    // parse the text into the root object
    return io::read_json(root, content);
}

/// simple error-reporting derserialize
bool io::deserialize_field(const Json::Value &val,
    json_serializable_t &object, std::string field)
{
    if (!object.deserialize(val[field]))
    {
        cerr << "Error, failed to deserialize (or missing) \"" + field + "\""
             << " field in input JSON." << endl;
        return false;
    }

    return true;
}

bool io::get_json_as_type(const Json::Value &obj,
    io::json_serializable_t &value)
{
    return value.deserialize(obj);
}

// JSON conversion overloads, note the terrible use of the comma operator
bool io::get_json_as_type(const Json::Value &obj, std::string &value) {
    return obj.isConvertibleTo(Json::stringValue) &&
        (value = obj.asString(), true);}

bool io::get_json_as_type(const Json::Value &obj, Json::Int &value) {
    return obj.isConvertibleTo(Json::intValue) &&
        (value = obj.asInt(), true);}

bool io::get_json_as_type(const Json::Value &obj, Json::UInt &value) {
    return obj.isConvertibleTo(Json::uintValue) &&
        (value = obj.asUInt(), true);}

bool io::get_json_as_type(const Json::Value &obj, Json::Int64 &value) {
    return obj.isConvertibleTo(Json::intValue) &&
        (value = obj.asInt64(), true);}

bool io::get_json_as_type(const Json::Value &obj, Json::UInt64 &value) {
    return obj.isConvertibleTo(Json::uintValue) &&
        (value = obj.asUInt64(), true);}

bool io::get_json_as_type(const Json::Value &obj, float &value) {
    return obj.isConvertibleTo(Json::realValue) &&
        (value = obj.asFloat(), true);}

bool io::get_json_as_type(const Json::Value &obj, double &value) {
    return obj.isConvertibleTo(Json::realValue) &&
        (value = obj.asDouble(), true);}

bool io::get_json_as_type(const Json::Value &obj, bool &value) {
    return obj.isConvertibleTo(Json::booleanValue) &&
        (value = obj.asBool(), true);}
