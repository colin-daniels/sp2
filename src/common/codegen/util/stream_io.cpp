#include "stream_io.hpp"
#include <limits>

void write_array(std::ostream &out, const std::vector<double> &data,
    const std::string prefix, const std::string suffix)
{
    auto old_precision = out.precision(
        std::numeric_limits<double>::max_digits10
    );

    out << "{" << suffix;
    for (std::size_t i = 0; i < data.size(); ++i)
    {
        out << prefix << data[i];
        if (i + 1 != data.size())
            out << ", ";

        out << suffix;
    }
    out << "}";

    out.precision(old_precision);
}
