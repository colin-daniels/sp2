#ifndef SP2_INCLUDE_JSON_JSON_H
#define SP2_INCLUDE_JSON_JSON_H
// Stand-in include file for json/json.h, which has had a history of
//  publishing code that triggers compiler warnings.

#include <boost/predef.h>

// disable warnings
#ifdef BOOST_COMP_GNUC
#pragma GCC system_header // note: effect extends to end of this file
#endif

#include <json/json.h>
#endif