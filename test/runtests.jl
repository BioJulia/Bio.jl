
function get_bio_fmt_specimens()
    path = Pkg.dir("Bio", "test", "BioFmtSpecimens")
    if !isdir(path)
        run(`git clone --depth 1 https://github.com/BioJulia/BioFmtSpecimens.git $(path)`)
    end
end

include("align/test_align.jl")
include("phylo/test_phylo.jl")
include("intervals/test_intervals.jl")
include("seq/test_seq.jl")
include("services/test_services.jl")
include("tools/test_tools.jl")
