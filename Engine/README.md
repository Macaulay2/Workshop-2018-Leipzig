# Engine

This is the collaboration area for the engine project at the 2018 Macaulay2
workshop in Leipzig.

## Team members

 * Dylan Peifer (djp282@cornell.edu)

## Potential projects

 1. Producing a set of documentation and examples for writing code in the
    engine and linking with a package. The goal is to make it possible for
    package authors to read about how and why (or why not) to move parts of
    their package into the engine. To start, a quick copy of some notes I took
    when writing a division algorithm function in the engine is in the
    `documentation` subdirectory.

 2. Writing engine functions to perform polynomial long division,
    minimalization, and interreduction. These functions should be heavily
    tested and benchmarked. The goal is to give top-level users more
    fine-grained control over the steps in producing a Groebner basis, and
    provide building blocks for a Groebner walk algorithm in the engine. A
    basic start to a package providing an interface to these functions is in
    the `reduction` subdirectory. You will not be able to use
    `rawDivisionAlgorithm` yet, but the other functions give simple top-level
    implementations of the operations we want.

 3. Writing an M2 kernel for Jupyter. Clark Spessard (not attending the
    workshop) has a [working version][1] that we could modify or improve. The
    goal is to distribute a working kernel for Jupyter with a standard
    Macaulay2 install.

Feel free to comment or suggest further projects!

[1]: https://github.com/clarkbsp/macaulay2-kernel

