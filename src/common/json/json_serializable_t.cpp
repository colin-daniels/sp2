#include "common/json/json_serializable_t.hpp"
#include "common/json/json.hpp"

using namespace sp2;

io::json_serializable_t::operator Json::Value() const
{
    Json::Value value;
    serialize(value);
    return value;
}

io::json_serializable_t &io::json_serializable_t::operator=(
    const Json::Value &obj)
{
    deserialize(obj);
    return *this;
}
