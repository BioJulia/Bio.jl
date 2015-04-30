
# Contributing

We welcome contributions in the form of pull requests. For your code to be
considered it must meet the following guidelines.

  * By making a pull request, you're agreeing to license your code under an MIT
    license. See LICENSE.md.

  * Types and functions must be documented using
    [Docile](https://github.com/MichaelHatherly/Docile.jl) style docstring.
    Documentation regarding specific implementation details that aren't relevent
    to users should be in the form of comments.

    Documentation may be omitted if the function not exported (i.e. only used
    internally) and is short and obvious. E.g. `cube(x) = x^3`.

  * All significant code must be tested. Tests should be organized into
    contexts, and into separate files based on module.

  * Contributions are included if the code has been reviewed by at least two
    team members and there is general consensus (or general lack of objections)
    that it's useful and fits with the intended scope of Bio.jl.

  * Code must be consistent with the prevailing style in Bio.jl, which includes,
    but is not necessarily limited to the following style guide.


## Style


  * Indent with 4 spaces.

  * Type names are camel case, with the first letter capitalized. E.g.
    `SomeVeryUsefulType`.

  * Function names, apart from constructors, are all lowercase. Include
    underscores between words only if the name would be hard to read without.
    E.g.  `start`, `stop`, `findletter` `find_last_digit`.

  * Generally try to keep lines below 80-columns, unless splitting a long line
    onto multiple lines makes it harder to read.

  * Files that declare modules should only declare the module, and import any
    modules that it requires. Any code should
    be included from separate files. E.g.

    ```julia
    module AwesomeFeatures

    using Compat, JSON

    include("feature1.jl")
    include("feature2.jl")

    end
    ```

  * Separate logical blocks of code with one blank line, or two blank lines for
    function/type definitions.

  * When extending method definitions, explicitly import the method.

    ```julia
    import Base: start, next, done
    ```

## Conduct

We adhere to the Julia [community standards](http://julialang.org/community/standards/).


