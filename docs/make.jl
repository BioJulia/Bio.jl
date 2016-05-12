using Documenter, Bio

custom_deps() = run(`pip install --user pygments mkdocs mkdocs-material`)
@osx? makedocs(doctest = false) : makedocs()
deploydocs(
           repo = "github.com/BioJulia/Bio.jl.git",
           julia = "0.4",
           osname = "linux",
           deps = custom_deps
           )
