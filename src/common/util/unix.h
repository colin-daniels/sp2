//
// Created by lampam on 7/12/17.
//

#ifndef SP2_UTIL_UNIX_H
#define SP2_UTIL_UNIX_H

#include <string>
#include <utility>

namespace sp2 {
namespace util {
namespace unix {

std::pair<std::string, std::string> split_path(std::string path);

}
}
}

#endif //SP2_UTIL_UNIX_H
