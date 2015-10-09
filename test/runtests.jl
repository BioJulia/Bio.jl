
function get_bio_fmt_specimens()
    path = Pkg.dir("Bio", "test", "BioFmtSpecimens")
    if !isdir(path)
        run(`git clone --depth 1 https://github.com/BioJulia/BioFmtSpecimens.git $(path)`)
    end
end

include("align/TestAlign.jl")
include("phylo/TestPhylo.jl")
include("intervals/TestIntervals.jl")
include("seq/TestSeq.jl")
include("services/TestServices.jl")
include("tools/TestTools.jl")
