//
// Created by lampam on 7/15/17.
//

// This is designated as the unique compilation unit which will define the
// global object for the numpy array API.

#define SP2_ACTUALLY_DO_PLEASE_IMPORT_ARRAY

#include "common/python/include_numpy.hpp"

int initialize_numpy()
{
    // import_array() is strange.
    //  - It is a macro with a return statement.
    //    But only on the failing branch!
    //  - Its expansion appears to trigger -Wconversion-null?
    //  - To make matters worse, its return type is either void or integral based on the python version.
    //
    // All of this absurd nonsense seems to be because it is intended to be used inside "your extension
    //  module's init function." But this, of course, does not apply to us.
    //
    // So instead we'll just call the undocumented function that the macro delegates to,
    // which is far more reasonably behaved.

    return _import_array();
}
