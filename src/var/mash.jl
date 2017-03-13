# MinHash
# ========
#
# Functions to get MASH distance for MinHash sketches
#
# see DOI: 10.1186/s13059-016-0997-x
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

function jaccarddistance(sketch1::MinHashSketch, sketch2::MinHashSketch)
    d = length(setdiff(sketch1.sketch, sketch2.sketch))
    l = length(sketch1)
    return (l-d) / (l+d)
end

function mashdistance(k::Int, j::Float64)
    return -1/k * log(2j / (1+j))
end

"""
    mashdistance(sketch1::MinHashSketch, sketch2::MinHashSketch)

Determines the MASH distance between two MinHash sketches. Requires that each
sketch is the same size, using kmers of the same length
"""
function mashdistance(sketch1::MinHashSketch, sketch2::MinHashSketch)
    sketch1.kmersize == sketch2.kmersize || error("sketches must have same kmer length")
    length(sketch1) == length(sketch2) || error("sketches must be the same size")

    j = jaccarddistance(sketch1, sketch2)
    k = sketch1.kmersize
    return mashdistance(k, j)
end
