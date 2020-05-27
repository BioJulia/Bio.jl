using Documenter, Bio, Pkg

makedocs(
    format = Documenter.HTML(),
    sitename = "Bio.jl",
    authors = replace(join(Pkg.TOML.parsefile("Project.toml")["authors"], ", "), r" <.*?>" => "" ),
    pages = [
        "Home"    => "index.md"
    ],
    
)

deploydocs(
    repo = "github.com/BioJulia/Bio.jl.git",
    push_preview = true,
    deps = nothing,
    make = nothing
)