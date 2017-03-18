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
    sketch1.kmersize == sketch2.kmersize || error("sketches must have same kmer length")
    length(sketch1) == length(sketch2) || error("sketches must be the same size")

    matches = 0
    sketchlen = length(sketch1)
    i = 1
    j = 1

    while i <= sketchlen && j <= sketchlen
        if sketch1.sketch[i] == sketch2.sketch[j]
            matches += 1
            i += 1
            j += 1
        elseif sketch1.sketch[i] < sketch2.sketch[j]
            while i <= sketchlen && sketch1.sketch[i] < sketch2.sketch[j]
                i += 1
            end
        elseif sketch2.sketch[j] < sketch1.sketch[i]
            while j <= sketchlen && sketch2.sketch[j] < sketch1.sketch[i]
                j += 1
            end
        end
    end

    if matches == sketchlen
        return 1.0
    else
        return matches / (2 * sketchlen - matches)
    end
end

function mashdistance(k::Integer, j::Float64)
    return -1/k * log(2j / (1+j))
end

"""
    mashdistance(sketch1::MinHashSketch, sketch2::MinHashSketch)

Determines the MASH distance between two MinHash sketches. Requires that each
sketch is the same size, using kmers of the same length
"""
function mashdistance(sketch1::MinHashSketch, sketch2::MinHashSketch)
    j = jaccarddistance(sketch1, sketch2)
    k = sketch1.kmersize
    return mashdistance(k, j)
end
