#ifndef SP2_JSON_SERIALIZABLE_HPP_HPP
#define SP2_JSON_SERIALIZABLE_HPP_HPP

#include <json/forwards.h>

namespace sp2 {
namespace io {

/// abstract base class for objects which can be serialized into/from json
class json_serializable_t
{
public:
    json_serializable_t() = default;
    json_serializable_t(const json_serializable_t&) = default;

    json_serializable_t& operator=(const Json::Value&);
    operator Json::Value() const;

    virtual ~json_serializable_t() = default;

    /// serialize object to Json::Value
    virtual bool serialize(Json::Value &output) const = 0;
    /// deserialize object from Json::Value
    virtual bool deserialize(const Json::Value &input) = 0;
};

} // namespace io
} // namespace sp2

#endif // SP2_JSON_SERIALIZABLE_HPP_HPP
