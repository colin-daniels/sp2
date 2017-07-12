//
// Created by lampam on 7/12/17.
//

#include "unix.h"

/*
#include "templates.hpp"

// FIXME holy moses this is a travesty of hacks upon hacks.
//       at this point we might as well use Boost::FileSystem

//---------------------------
// HACK forward-define delegated stumps before importing libgen
//      because libgen may define basename or dirname as macros,
//      preventing us from defining our own correspondingly-named functions.
static std::string basename_impl(const std::string &path);
static std::string dirname_impl(const std::string &path);

std::string sp2::util::unix::basename(const std::string &path) { return basename_impl(path); }
std::string sp2::util::unix::dirname(const std::string &path)  { return dirname_impl(path); }
//---------------------------

#include "libgen.h"
 */

using namespace std;

static string& lstrip(string &s, char c) {
    size_t i = s.find_first_not_of(c);
    if (i != string::npos) {
        s.erase(0, i);
    }
    return s;
}

static string& rstrip(string &s, char c) {
    while (s.size() > 0 && s.rfind(c) == s.size() - 1) {
        s.pop_back();
    }
    return s;
}

pair<string, string> sp2::util::unix::split_path(string path) {
    // ""
    if (path.empty()) { return make_pair(".", "."); }

    // strip trailing /
    rstrip(path, '/');

    // "/", "//", "///", ...
    if (path.empty()) { return make_pair("/", "."); }

    size_t last_slash = path.rfind("/");

    // "just-a-filename"
    if (last_slash == string::npos) { return make_pair(".", path); }

    // "/filename"
    if (last_slash == 0) { return make_pair("/", path.erase(0,1)); }

    // General case. Examples:
    //   FUNCTION INPUT              DIRNAME      BASENAME
    //   some/dir/basename       ->  some/dir     basename
    //   /some/dir/basename      ->  /some/dir    basename
    //   some/dir/basename/      ->  some/dir     basename
    //   some/dir/.              ->  some/dir     .
    //   ///silly/dir/basename   ->  ///silly/dir basename
    //   /silly/dir///basename   ->  /silly/dir   basename
    auto mid = path.begin() + last_slash;
    auto dirname = string(path.begin(), mid);
    rstrip(dirname, '/');
    auto basename = string(mid+1, path.end());
    return make_pair(dirname, basename);
};

/*
// FIXME not thread safe due to how the libgen function may store its result?
string basename_impl(const string &path) {
    // some notes:
    //  * we need a mutable char * because the libgen function may trash its input string
    //  * we need to copy the output ASAP as it might point to static memory which will only temporarily be valid.
    char *in = new char[path.size()];
    auto guard = sp2::scope_guard([&] { delete [] in; });
    copy(path.begin(), path.end(), in);

    return string(basename(in));
}

// FIXME not thread safe due to how the libgen function may store its result?
string dirname_impl(const string &path) {
    // some notes:
    //  * we need a mutable char * because the libgen function may trash its input string
    //  * we need to copy the output ASAP as it might point to static memory which will only temporarily be valid.
    char *in = new char[path.size()];
    auto guard = sp2::scope_guard([&] { delete [] in; });
    copy(path.begin(), path.end(), in);

    return string(dirname(in));
}
*/