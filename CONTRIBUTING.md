
# Contributing

We welcome contributions in the form of pull requests. For your code to be
considered it must meet the following guidelines.

  * By making a pull request, you're agreeing to license your code under an MIT
    license. See LICENSE.md.

  * Types and functions must be documented using Julia's [docstrings](http://docs.julialang.org/en/latest/manual/documentation/).
    Documentation regarding specific implementation details that aren't relevant
    to users should be in the form of comments.

    Documentation may be omitted if the function is not exported (i.e. only used
    internally) and is short and obvious. E.g. `cube(x) = x^3`.

  * All significant code must be tested. Tests should be organized into
    contexts, and into separate files based on module.

  * Contributions are included if the code has been reviewed by at least two
    team members who are **not** the author of the proposed contribution,
    and there is general consensus (or general lack of objections) that it's useful
    and fits with the intended scope of Bio.jl.

  * Code must be consistent with the prevailing style in Bio.jl, which includes,
    but is not necessarily limited to the following style guide.

  * Code contributed should be compatible with Julia v0.4.


## Style


  * Indent with 4 spaces.

  * Type names are camel case, with the first letter capitalized. E.g.
    `SomeVeryUsefulType`.

  * Module names are also camel case.

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

    using IntervalTrees, JSON

    include("feature1.jl")
    include("feature2.jl")

    end
    ```
  * Files that declare modules should have the same name name of the module, e.g
    the module `SomeModule` is declared under the file `SomeModule.jl`.

  * Separate logical blocks of code with one blank line, or two blank lines for
    function/type definitions.

  * When extending method definitions, explicitly import the method.

    ```julia
    import Base: start, next, done
    ```

  * Document functions using bare docstrings before a definition:

  ```julia
  "This function foo's something"
  foo(x) = 2*x
  ```

  * Functions that get or set variables in a type should not be prefixed with 'get' or 'set'. The getter should be named for the variable it sets, and the setter should have the same name as the getter, with the suffix `!`. For exmaple, for the variable `names`:

  ```julia
  name(node) # get node name
  name!(node, "somename") # set node name
  ```

## Conduct

We adhere to the Julia [community standards](http://julialang.org/community/standards/).
