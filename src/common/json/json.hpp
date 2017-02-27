//
// Created by cc on 9/7/16.
//

#ifndef SP2_JSON_HPP
#define SP2_JSON_HPP

#include "common/json/json_serializable_t.hpp"
#include "common/enums.hpp"

#include <string>
#include <utility>
#include <iostream>
#include <vector>

#include <json/json.h>

namespace sp2 {
namespace io {

////////////////////////////////////////////////////////////////////////////////
// JSON IO Functions.                                                         //
////////////////////////////////////////////////////////////////////////////////

/// write json value to string
bool write_json(const Json::Value &root, std::string &output,
    bool styled = true);
/// write json to a file
bool write_json_file(const Json::Value &root, std::string filename,
    bool styled = true);

/// read json value from string
bool read_json(Json::Value &root, const std::string &input);
/// read json from a file
bool read_json_file(Json::Value &root, std::string filename);

/// simple error-reporting derserialize
bool deserialize_field(const Json::Value &val,
    json_serializable_t &object, std::string field);

////////////////////////////////////////////////////////////////////////////////
// Utilities                                                                  //
////////////////////////////////////////////////////////////////////////////////

bool get_json_as_type(const Json::Value &obj, json_serializable_t &value);
bool get_json_as_type(const Json::Value &obj, std::string &value);
bool get_json_as_type(const Json::Value &obj, Json::Int &value);
bool get_json_as_type(const Json::Value &obj, Json::UInt &value);
bool get_json_as_type(const Json::Value &obj, Json::Int64 &value);
bool get_json_as_type(const Json::Value &obj, Json::UInt64 &value);
bool get_json_as_type(const Json::Value &obj, float &value);
bool get_json_as_type(const Json::Value &obj, double &value);
bool get_json_as_type(const Json::Value &obj, bool &value);

template<
    typename T,
    typename = std::enable_if_t<
        std::is_convertible<const Json::Value&, T>::value
    >
>
bool get_json_as_type(const Json::Value &obj, T& value)
{
    value = obj;
    return true;
}

template<
    typename T,
    typename = std::enable_if_t<
        std::is_enum<T>::value
    >,
    std::size_t N = 2
>
bool get_json_as_type(const Json::Value &obj, T& value)
{
    if (!obj.isConvertibleTo(Json::stringValue))
        return false;

    value = enum_from_str<std::remove_cv_t<T>>(obj.asString());
    return true;
}

template<
    typename T,
    std::size_t N
>
bool get_json_as_type(const Json::Value &obj, T (&arr)[N])
{
    if (!obj.isArray() || obj.size() != N)
        return false;

    bool ok = true;
    for (int i = 0; i < N; ++i)
        if (!get_json_as_type(obj[i], arr[i]))
            ok = false;

    return ok;
}

template<
    typename T,
    std::size_t N,
    std::size_t M
>
bool get_json_as_type(const Json::Value &obj, T (&arr)[N][M])
{
    if (!obj.isArray() || obj.size() != N)
        return false;

    bool ok = true;
    for (int i = 0; i < N; ++i)
        if (!get_json_as_type(obj[i], arr[i]))
            ok = false;

    return ok;
}

template<typename T>
bool get_json_as_type(const Json::Value &obj, std::vector<T> &vec)
{
    if (!obj.isArray())
        return false;

    vec.clear();

    bool ok = true;
    for (int i = 0; i < obj.size(); ++i)
    {
        vec.emplace_back();
        if (!get_json_as_type(obj[i], vec.back()))
            ok = false;
    }

    return ok;
}


template<
    typename T,
    typename = std::enable_if_t<
        std::is_convertible<const T&, Json::Value>::value
    >
>
void get_type_as_json(const T& value, Json::Value &obj)
{
    obj = value;
}


template<
    typename T,
    typename = std::enable_if_t<
        std::is_enum<T>::value
    >,
    std::size_t N = 2
>
void get_type_as_json(const T& value, Json::Value &obj)
{
    obj = enum_to_str(value);
}


template<typename T>
void get_type_as_json(const std::vector<T> &value, Json::Value &obj)
{
    // create empty array first
    obj = Json::Value(Json::arrayValue);

    for (auto &elem : value)
    {
        Json::Value temp;
        get_type_as_json(elem, temp);
        obj.append(temp);
    }
}

template<
    typename T,
    std::size_t N
>
void get_type_as_json(const T (&arr)[N], Json::Value &obj)
{
    obj = Json::Value(Json::arrayValue);
    for (std::size_t i = 0; i < N; ++i)
    {
        Json::Value temp;
        get_type_as_json(arr[i], temp);
        obj.append(temp);
    }
}

template<
    typename T,
    std::size_t N,
    std::size_t M
>
void get_type_as_json(const T (&arr)[N][M], Json::Value &obj)
{
    obj = Json::Value(Json::arrayValue);
    for (std::size_t i = 0; i < N; ++i)
    {
        Json::Value temp;
        get_type_as_json(arr[i], temp);
        obj.append(temp);
    }
}

////////////////////////////////////////////////////////////////////////////////
// Templates.                                                                 //
////////////////////////////////////////////////////////////////////////////////

////////////////////
// serialization
inline void serialize_basic(Json::Value &output) {}

template<typename T, typename ...Args>
void serialize_basic(Json::Value &output,
    const std::string &key, T&& value,
    Args &&...args_left)
{
    get_type_as_json(std::forward<T>(value), output[key]);
    serialize_basic(output, std::forward<Args>(args_left)...);
}

////////////////////
// deserialization
inline bool deserialize_basic(const Json::Value &input) {return true;}

template<typename T, typename ...Args>
bool deserialize_basic(const Json::Value &input,
    const std::string &key, T&& value, Args &&...args_left)
{
    // TODO: Output warning when failing to deserialize a present field
    bool a = input.isMember(key) &&
             get_json_as_type(input[key], std::forward<T>(value));

    bool b = deserialize_basic(input, std::forward<Args>(args_left)...);

    return a && b;
}

/// \brief read an enumeration-like field
/// \tparam T enumeration-like type determined from the field
/// \param value Json::Value& input JSON object
/// \param field_name the field name to be checked in the JSON input
/// \param choices std::vector<std::pair<std::string, T>> list of
///        acceptable inputs, maps from string to value
/// \param default_choice T the value to be returned if no match is
///        found or if the function fails to read the field
/// \param output bool whether to output to standard output on error
template<class T>
T deserialize_enum(const Json::Value &value, std::string field_name,
    std::vector<std::pair<std::string, T>> choices,
    T default_choice, bool output = false)
{
    if (choices.empty())
        return default_choice;

    std::string enum_string = value.get(field_name, "").asString();
    if (enum_string.empty())
    {
        if (output)
            std::cout << "Error, failed to deserialize (or missing) "
                      << "required " << field_name << " in input JSON."
                      << std::endl;
        return default_choice;
    }

    // ignore case for conversion from string->enum
    for (char &c : enum_string)
        c = static_cast<char>(::tolower(c));

    // check if the string matches any input names
    for (auto &choice : choices)
        if (enum_string == choice.first)
            return choice.second;

    // if no match was found
    if (output)
    {
        std::cout << "Error, incorrect (or no) " << field_name
                  << " specified in input JSON. Accepted inputs are: "
                  << '\"' << choices[0].first << '\"';
        for (std::size_t i = 1; i < choices.size(); ++i)
            std::cout << ", \"" << choices[i].first << '\"';

        std::cout << std::endl;
    }
    return default_choice;
}

} // namespace io
} // namespace sp2

#endif //SP2_JSON_HPP
