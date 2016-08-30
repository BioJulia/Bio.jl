# BAM
# ===
#
# The BAM file format.
#
# http://samtools.github.io/hts-specs/SAMv1.pdf

"""
The BAM file format.
"""
immutable BAM <: Bio.IO.FileFormat end

include("bai.jl")
include("record.jl")
include("reader.jl")
include("writer.jl")
include("auxdict.jl")
include("intersect.jl")

function Base.open(filename::AbstractString, mode::AbstractString, ::Type{BAM})
    if mode[1] == 'r'
        return open(filename, BAM)
    elseif mode[1] == 'w'
        return BAMWriter(BGZFStream(filename, mode))
    end
    error("invalid open mode")
end
