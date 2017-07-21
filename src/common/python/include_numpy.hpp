//
// Created by lampam on 7/15/17.
//

//--------------------------
// TL;DR:
//
//  * Every compilation unit in 'sp2' which wants to import numpy headers
//    and which is not part of an extension module should import this instead.
//
//    (extension modules should instead follow the recommendations here:
//     https://docs.scipy.org/doc/numpy-1.10.1/reference/c-api.array.html#miscellaneous)
//
//  * initialize_numpy should be called early during program execution.
//    If it isn't, you'll get segfaults on innocent looking numpy API calls.
//
//--------------------------

// Importing numpy.  Oh, this is fun.

// numpy has a 'static' ARRAY_API object which must be initialized using a function
// called import_array.  You can read about it here, but their recommendations must be
// taken with a grain of salt since they assume you're writing an extension module:
//
// https://docs.scipy.org/doc/numpy-1.10.1/reference/c-api.array.html#miscellaneous

// Defining PY_ARRAY_UNIQUE_SYMBOL replaces the static object with a non-static global,
//  allowing multiple compilation units to share the same API object...

#define  PY_ARRAY_UNIQUE_SYMBOL  sp2_common_ARRAY_API

// ...but now we must contend with the One Definition Rule.
//
// To help cope with this, numpy provides NO_IMPORT_ARRAY to disable the API object
//  definition, which should be specified in every compilation unit except one.
//
// That, of course, is an absurd default, so we invert it.

#ifndef SP2_ACTUALLY_DO_PLEASE_IMPORT_ARRAY
#define NO_IMPORT_ARRAY
#endif

// Prohibit the use of APIs deprecated in numpy 1.10.
// We will force all consumers of this header to share this, because it appears
//  that this may actually impact the ABI (due to the way they hid deprecated data members).
// I am not sure.

#define NPY_NO_DEPRECATED_API NPY_1_10_API_VERSION

// Set up numpy. Call once during program invocation.
// If you forget it, expect segfaults.
//
// Returns a negative value and sets a Python exception (PyErr API) on failure.
int initialize_numpy();

#include "numpy/arrayobject.h"
