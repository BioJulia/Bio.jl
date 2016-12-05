using Documenter, Bio

makedocs()
deploydocs(
    deps = Deps.pip("mkdocs", "pygments", "mkdocs-biojulia"),
    repo = "github.com/BioJulia/Bio.jl.git",
    julia = "0.5",
    osname = "linux",
)
