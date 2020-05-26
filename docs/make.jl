using Documenter, Bio, Pkg

makedocs(
    format = Documenter.HTML(),
    sitename = "Bio.jl",
    authors = replace(join(Pkg.TOML.parsefile("Project.toml")["authors"], ", "), r" <.*?>" => "" ),
    pages = [
        "Home"    => "index.md",
        "Install" => "man/install.md",
        "Bio.Seq" => "man/seq.md",
        "Bio.Align" => "man/alignments.md",
        "Bio.Intervals" => "man/intervals.md",
        "Bio.Var" => "man/var.md",
        "Bio.Structure" => "man/structure.md",
        "Bio.Util" => "man/util.md",
        "Bio.Tools" => "man/tools.md",
        "IO API" => "man/reading.md"
    ],
    
)

deploydocs(
    repo = "github.com/BioJulia/Bio.jl.git",
    push_preview = true,
    deps = nothing,
    make = nothing
)