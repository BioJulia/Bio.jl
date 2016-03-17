using Lapidary, Bio

custom_deps() = run(`pip install --user pygments mkdocs mkdocs-material`)
makedocs()
cp(
   Pkg.dir("Bio", "docs/src/assets/images/biojl.png"),
   Pkg.dir("Bio", "docs/build/assets/images/biojl.png")
   )
deploydocs(
           repo = "github.com/BioJulia/Bio.jl.git",
           julia = "0.4",
           osname = "linux",
           deps = custom_deps
           )
