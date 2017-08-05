
`to_python` and `from_python`
=============================

`to_python` turns a variety of C++ data types into python equivalents.

`from_python` reads python data back into strongly typed C++.

`to_python`
-----------

```c++
// Approximate signature:

template<typename C> bool to_python(const py_scoped_t& py, C &c);
template<typename C> bool to_python(const py_object_t& py, C &c);
```

`to_python` generally tries to produce types which approximate or represents
the C++ type in spirit.  For cases that are unclear, selection of the
desired python type may be accomplished through the use of wrapper types.
(this is conjecture; these wrappers don't exist yet)

`to_python` is generally unlikely to fail, but when it does, all references
created by it are guaranteed to be cleaned up, and it will leave the output
argument NULL and return `false`.

`from_python`
-------------

```c++
// Approximate signature:

template<typename C> bool from_python(const C& c, py_scoped_t &py);
template<typename C> bool from_python(const C& c, py_object_t &py);
```

When both are implemented, `from_python` must be able to recover an object
from the output of `to_python`, but there is no general restriction beyond that.
`from_python` implementations tend to be fairly permissive in what they accept,
performing many of the sorts of conversions accepted in common Python practice
(e.g. accepting an int or Fraction in place of a float, or an iterable in place
of a list, etc.).

On failure, `from_python` returns `false`, and leaves the output argument in
an unspecified state.

Failure is not uncommon, and in our case generally results from poor output
returned by a user script; hence good error messages are a must.
The current implementations typically try to take advantage of python's
own fantastic error messages and tracebacks by using `PyErr_Print`
before returning false.

Of course, unconditionally printing errors makes the current implementations
unsuitable for use cases such as "attempt to deserialize as type A, then
attempt to deserialize as type B if that fails". This may be addressed in
the future.

`py_scoped_t` versus `py_object_t`
----------------------------------

All conversions are implemented for both `py_scoped_t` and `py_object_t`.
Use whatever type is available to the code you are writing.

`py_object_t` has convenience functions built into it that don't take output
arguments, and instead throw exceptions. (Er, that's, uh... in addition to
unconditionally printing to stderr. Sorry!)

Code internal to the python subdirectory will probably tend to use the
py_scoped_t versions, since these are (probably) significantly faster for
generic container types.