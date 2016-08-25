using Documenter, Bio

custom_deps() = run(`pip install --user pygments mkdocs mkdocs-biojulia`)
makedocs()
deploydocs(
           repo = "github.com/BioJulia/Bio.jl.git",
           julia = "0.5",
           osname = "linux",
           deps = custom_deps
           )
